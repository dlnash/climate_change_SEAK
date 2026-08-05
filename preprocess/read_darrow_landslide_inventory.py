import pandas as pd

# Read the file
df = pd.read_excel("../out/darrow_AK_landslide_inventory.xlsx", header=None)

# Rename columns from the two header rows
df.columns = [
    "Date",
    "Location",
    "Latitude",
    "Region",
    "Cause",
    "Injuries_Fatalities",
    "Damage_Impact",
    "Source",
]

# Remove the two header rows
df = df.iloc[2:].reset_index(drop=True)

records = []

# Every record occupies two rows
for i in range(0, len(df), 2):

    row1 = df.iloc[i]
    row2 = df.iloc[i + 1]
    records.append({
            "Date": pd.to_datetime(str(row1["Date"]).rstrip(","), errors="coerce"),
            "Julian_DOY": pd.to_numeric(row2["Date"], errors="coerce"),
            "Location": row1["Location"],
            "Latitude": pd.to_numeric(str(row1["Latitude"]).rstrip(","), errors="coerce"),
            "Longitude": pd.to_numeric(row2["Latitude"], errors="coerce"),
            "Region": row1["Region"],
            "Cause": row1["Cause"],
            "Injuries_Fatalities": row1["Injuries_Fatalities"],
            "Damage_Impact": row1["Damage_Impact"],
            "Source": row1["Source"],
})

new_df = pd.DataFrame(records)

new_df["Date"] = pd.to_datetime(new_df["Date"])

# Define filter criteria
start_date = pd.Timestamp("1981-01-01")
end_date   = pd.Timestamp("2019-12-31")

filtered_df = new_df[
    # Date range
    (new_df["Date"] >= start_date) &
    (new_df["Date"] <= end_date) &

    # Cause
    (new_df["Cause"].isin(["R", "R*", "F/T"])) &

    # Latitude (54°N to 60°N)
    (new_df["Latitude"].between(54, 60)) &

    # Longitude (130°W to 140°W)
    # Note: western longitudes are negative
    (new_df["Longitude"].between(-140, -130))
].copy()

print(f"Number of events: {len(filtered_df)}")
filtered_df