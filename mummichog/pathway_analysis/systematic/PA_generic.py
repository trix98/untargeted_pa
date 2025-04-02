# For general data science and matrix manipulation
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# For pathway analysis in python
import sspa

# For plotting
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px

def load_mummichog_output_data(mummichog_output_filename):
    '''
    filename e.g., "mummichog/runs/trans_omic_covid_data.rsd_1_default_p/tables/userInput_to_EmpiricalCompounds.tsv"
    '''

    mum_out_map = pd.read_csv(mummichog_output_filename, sep='\t')
    print(f"Of the {len(mum_out_map)} features/compounds (e.g. M70.065...) in mummichog output {len(np.unique(mum_out_map['CompoundID_from_user']))} are unique.")

    drop_cols = [col for col in mum_out_map.columns if col not in ["compound_names", "CompoundID_from_user"]]
    mum_out_map.drop(columns=drop_cols, inplace=True)
    mum_out_map = mum_out_map.rename(columns={'CompoundID_from_user': 'namecustom'})
    mum_out_map['annotation'] = [str(row['compound_names']).split(';')[0] for _, row in mum_out_map.iterrows()]
    mum_out_map.drop(columns=['compound_names'], inplace=True)

    return mum_out_map

def merge_meta(mum_output, mum_input_filename, meta_filename):
    
    mum_input = pd.read_csv(mum_input_filename)
    drop_cols = [col for col in mum_input.columns if col[0]!="2" and col not in ["namecustom"]]
    mum_input.drop(columns=drop_cols, inplace=True)
    # print(mum_output.head())
    # print(len(mum_output), len(mum_output.columns))
    # print("------------")
    # print(mum_input.head())
    # print(len(mum_input), len(mum_input.columns))

    # print(len([el for el in mum_output['namecustom'] if el in mum_input['namecustom'].tolist()]))
    # print("------------")

    df_merged = mum_output.merge(mum_input, on='namecustom')

    df_merged = df_merged.drop(columns=['namecustom'])
    df_merged = df_merged.rename(columns={df_merged.columns[0]: 'sample_id'})
    df_merged = df_merged.set_index('sample_id')
    df_merged = df_merged.transpose()

    # print(len(df_merged), 'samples')
    # print(len(df_merged.columns), 'compounds')
    # print(df_merged.head())

    df_meta = pd.read_csv(meta_filename)
    df_meta = df_meta[['sample_id', 'is_bad', 'age']]
    df_meta = df_meta[df_meta['sample_id'].isin(df_merged.index)]
    df_meta = df_meta.set_index('sample_id')

    df_combined = df_merged.merge(df_meta, left_index=True, right_index=True)
    df_values = df_combined.iloc[:, :-2]
    df_meta = df_combined.iloc[:, -2:]

    return df_values, df_meta

def process_data(df_values):
    data_filt = df_values.loc[:, df_values.isin([' ', np.nan, 0]).mean() < 0.5]
    imputed_mat = data_filt.fillna(data_filt.median())
    log2_mat = np.log2(imputed_mat)
    processed_data = pd.DataFrame(StandardScaler().fit_transform(log2_mat), columns=imputed_mat.columns, index=imputed_mat.index)
    return processed_data

def do_pca(processed_data, df_meta):
    PCA_covid = PCA(n_components=2)
    PCA_scores = pd.DataFrame(PCA_covid.fit_transform(processed_data), columns=['PC1', 'PC2'], index=df_meta.index)
    PCA_scores['is_bad'] = df_meta['is_bad'].values
    PCA_scores['age'] = df_meta['age'].values

    sns.set_style('darkgrid')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    sns.scatterplot(
        data=PCA_scores,
        x='PC1',
        y='PC2',
        hue='is_bad',
        linewidth=0.5,
        edgecolor='k',
        s=50,
        # palette='Reds',
        # hue_order=['0', '1-2', '3-4', '5-7'],
        ax=ax1
    )

    sns.scatterplot(
        data=PCA_scores,
        x='PC1',
        y='PC2',
        hue='age',
        linewidth=0.5,
        edgecolor='k',
        s=50,
        palette='Reds',
        hue_order=['0', '1-2', '3-4', '5-7'],
        ax=ax2
    )
    plt.show()

def cpd_name_to_KEGG(processed_data, conversion_file):
    # processed_data.to_csv("transomic_covid_data_processed_gt.csv")
    compound_names = processed_data.columns.tolist()
    # df = pd.DataFrame({"Compounds": compound_names})
    # df.to_csv(f"Compounds_to_metaboanalyse_rsd_{suffix}.csv")

    conversion_table = pd.read_csv(conversion_file)

    print(len(compound_names), "compounds in the dataset")
    kegg_matches = conversion_table[(conversion_table["Comment"] == 1) & (conversion_table["KEGG"].isnull()==False)]["KEGG"]
    print(len(kegg_matches), "compounds from the dataset that have KEGG IDs")

    kegg_human_pathways  = sspa.process_kegg(organism="hsa")
    processed_data_mapped = sspa.map_identifiers(conversion_table, output_id_type="KEGG", matrix=processed_data)
    all_kegg_cpds = set(sum(sspa.utils.pathwaydf_to_dict(kegg_human_pathways).values(), []))
    print(f"All KEGG cpds: {len(all_kegg_cpds)}")
    mapped_annotated_cpds = set(processed_data_mapped.columns) & all_kegg_cpds
    print(len(mapped_annotated_cpds), "compounds present in both the dataset and kegg pathways")

    pathways_dict = sspa.utils.pathwaydf_to_dict(kegg_human_pathways)

    # How many pathways contain at least two mapped compounds?
    pathways_present = {k: v for k, v in pathways_dict.items() if len([i for i in processed_data_mapped.columns if i in v]) > 1}
    print(len(pathways_present), "pathways contain at least 2 mapped compounds")

    def plot_box():
        sns.set_context('notebook')
        sns.set_style('ticks')
        sns.barplot(
            y=[len(compound_names), len(kegg_matches), len(mapped_annotated_cpds)],
            x=['Original', 'Mapping to KEGG ID', 'Mapping to KEGG ID \n and annotated to KEGG']
            )
        plt.tight_layout()
        plt.show()
    
    def plot_funnel():
        data = dict(count=[len(compound_names), len(kegg_matches), len(mapped_annotated_cpds)],
            label=['Mummichog Output Annotations', 'Annotations with KEGG ID', 'Annotated to KEGG pathways'])

        fig = px.funnel(data, x='count', y='label')
        # fig.show(renderer="browser")
        fig.write_image("/Users/pranathipoojary/Projects/mummichog_proj/untargeted_pa/mummichog/pathway_analysis/systematic/rsd_1_funnel.png")
        
    plot_funnel()

    return processed_data_mapped, kegg_human_pathways

def ora_analysis(processed_data_mapped, kegg_human_pathways, df_meta, save_name):
    ora = sspa.sspa_ora(
    mat=processed_data_mapped, # Processed data matrix
    metadata=df_meta['is_bad'], # metadata column
    pathways=kegg_human_pathways, # pathway dataframe
    DA_cutoff=0.01, # t-test cutoff to select differential metabolites
    custom_background=None) # None sets to the default background set which are all annotated compounds provided in the input matrix

    # perform ORA
    ora_res = ora.over_representation_analysis()
    print('There are', len(ora.DA_molecules), 'differential metabolites')

    top_20_pathways = ora_res.sort_values(by="P-value").iloc[0:20, :]
    top_20_pathways["Pathway_name"] = top_20_pathways["Pathway_name"].str.replace(" - Homo sapiens \\(human\\)", "", regex=True)

    plt.figure(figsize=(8, 6))
    sns.set_style('ticks')
    sns.barplot(
        data=top_20_pathways,
        y="Pathway_name",
        x="P-value",
        orient="h",
        palette="viridis"
        )
    plt.axvline(0.05, c="black")
    plt.title('Top 20 KEGG Human pathways (ORA)')
    plt.xlim([0, 0.5])
    plt.tight_layout()
    plt.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.show()


def mummichog_pathways(csv_file, save_name="kegg_pathways.png"):
    df = pd.read_excel(csv_file)
    df.columns = df.columns.str.strip().str.lower()  # Convert to lowercase and remove spaces
    df.rename(columns={'p-value': 'P-value', 'pathway': 'Pathway_name'}, inplace=True)
    df = df.sort_values(by="P-value").iloc[0:20, :]

    plt.figure(figsize=(10, 8))
    sns.set_style('ticks')
    sns.barplot(
        data=df,
        y="Pathway_name",
        x="P-value",
        orient="h",
        palette="viridis"
    )
    plt.axvline(0.05, c="black", linestyle="--", label="Significance Threshold (0.05)")
    plt.title("KEGG Human Pathways (Mummichog's Pathway Predictions)")
    # plt.xlabel("P-value")
    plt.xlim([0, 0.5])
    # plt.ylabel("Pathway")
    # plt.legend()
    plt.tight_layout()
    
    plt.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.show()


def run_PA(mum_ouput_file, mum_input_file, meta_file, conversion_file, mummichog_PA, save_dir):
    
    
    mum_output = load_mummichog_output_data(mum_ouput_file) # 2 cols: features (M70.06...) and annotation (1-pyrroline)
    df_values, df_meta = merge_meta(mum_output, mum_input_file, meta_file)

    processed_data = process_data(df_values)

    # do_pca(processed_data, df_meta)

    processed_data_mapped, kegg_human_pathways = cpd_name_to_KEGG(processed_data, conversion_file)
    ora_analysis(processed_data_mapped, kegg_human_pathways, df_meta, f"{save_dir}ORA_top_20_rsd_1.png")
    mummichog_pathways(mummichog_PA, save_name=f"{save_dir}mummichog_PA_top_20_rsd_1.png")

if __name__ == "__main__":
    mum_ouput_file = "/Users/pranathipoojary/Projects/mummichog_proj/untargeted_pa/mummichog/pathway_analysis/systematic/mummichog_outputs/rsd_1_userInput_to_EmpiricalCompounds.tsv"
    mum_input_file = "/Users/pranathipoojary/Projects/mummichog_proj/untargeted_pa/mummichog/pathway_analysis/systematic/prs_feat_table_pos.csv"
    meta_file = "/Users/pranathipoojary/Projects/mummichog_proj/untargeted_pa/mummichog/pathway_analysis/systematic/meta_data.csv"
    conversion_file = "/Users/pranathipoojary/Projects/mummichog_proj/untargeted_pa/mummichog/pathway_analysis/systematic/name_maps/name_map_rsd_1.csv"
    save_dir = "/Users/pranathipoojary/Projects/mummichog_proj/untargeted_pa/mummichog/pathway_analysis/systematic/PA_results/"
    mummichog_PA = "/Users/pranathipoojary/Projects/mummichog_proj/untargeted_pa/mummichog/pathway_analysis/systematic/mummichog_pathways/rsd_1_mcg_pathwayanalysis.xlsx"
    run_PA(mum_ouput_file, mum_input_file, meta_file, conversion_file, mummichog_PA, save_dir)
