import os
import re

# ===== USER CONFIGURATION =====
# Path to the directory containing the files
directory_path = "/share/neutrino/snoplus/Data/FullFill_2p2/ratds/"

# Path to save the output text file
output_dir = "/lstore/sno/joankl/solar_analysis/real_data/"
output_fname = "run_numbers.txt"
full_output_dir = output_dir + output_fname

# ===== MAIN SCRIPT =====

# List to store the extracted run numbers
run_numbers = []

# === Regex pattern to match both formats:
# Matches:
#   Analysis20_r0000305600_s000_p000.root
#   Analysis20R_r0000308949_s009_p000.root
# Explanation:
#   ^Analysis20(R)?_r  → Matches "Analysis20_r" or "Analysis20R_r"
#   0*  → Matches any leading zeros
#   (\d+) → Captures the actual run number without leading zeros
pattern = re.compile(r"^Analysis20(R)?_r0*(\d+)_s\d+_p\d+\.root$")

# Loop through all files in the given directory
for filename in os.listdir(directory_path):
    # Check if the filename matches our pattern
    match = pattern.match(filename)
    if match:
        # Extract the run number (from the first capture group)
        run_number = match.group(2)
        run_numbers.append(run_number)

# ===== OUTPUT =====
# Print the list of run numbers
print("Extracted run numbers:", run_numbers)

# Save the run numbers into a .txt file (one per line)
with open(full_output_dir, "w") as f:
    for rn in run_numbers:
        f.write(rn + "\n")

print(f"Run numbers saved to: {full_output_dir}")
