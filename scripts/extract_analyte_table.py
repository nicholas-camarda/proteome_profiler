import os

import pandas as pd
import tabula as tb

# Define the pages to extract
pages_to_extract = [17, 18, 19]
filename = "cytoXL array kit - protocol.pdf"
# pages_to_extract = [15, 16]
# filename = "Mouse Angiogenesis Array Kit - Protocol.pdf"

# Remove the extension of the filename
filename_without_ext = os.path.splitext(filename)[0]
# Extract the tables from the PDF
tables = tb.read_pdf(filename, pages=pages_to_extract, multiple_tables=True)

# Combine the extracted tables into a single DataFrame
df = pd.concat(tables, ignore_index=True)

# Save the DataFrame to an Excel file
df.to_excel(f"output/{filename_without_ext}.xlsx", index=False)
