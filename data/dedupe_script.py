import pandas as pd

# 1. Define the input and output filenames
INPUT_FILE = "antigen_cdr3_pairs.csv"
OUTPUT_FILE = "antigen_cdr3_pairs_deduped.csv"

# 2. Load the data
# Make sure the file is in the same directory, or specify the full path.
df = pd.read_csv(INPUT_FILE)

# 3. Drop duplicates based on the primary sequence columns
# This keeps only one instance of each unique Antigen-CDR3 pair.
df_deduped = df.drop_duplicates(subset=['antigen_seq', 'cdr3_seq'], keep='first')

# 4. Save the new, clean file
df_deduped.to_csv(OUTPUT_FILE, index=False)

print(f"Original pairs: {len(df)}")
print(f"Unique pairs: {len(df_deduped)}")
print(f"Deduped file saved as: {OUTPUT_FILE}")