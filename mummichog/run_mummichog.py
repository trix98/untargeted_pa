import mummichog
import subprocess
import os

curr_dir = os.getcwd()
input_file = os.path.join(curr_dir, "mummichog/input_data/mummichog_input_ttest.txt")
work_dir = os.path.join(curr_dir, "mummichog/runs")
output_dir = "run_2_p_val"
cut_off_p = "0.05"

command = [
    "mummichog",
    "-f", input_file,
    "-o", output_dir,
    "--workdir", work_dir,
    "-m", "pos_default",
    "-c", cut_off_p,
]

try:
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    
    print("Mummichog execution completed successfully.")
    print("\n--- Terminal Output (stdout) ---\n")
    print(result.stdout) 
    print("\n--- Terminal Errors/Warnings (stderr) ---\n")
    print(result.stderr)

except subprocess.CalledProcessError as e:
    print("Error occurred while running mummichog:")
    print("\n--- Error Output (stderr) ---\n")
    print(e.stderr)
