import pandas as pd

#enriched
df["odds_ratio"] = (df["k"]/df["N"]) / (df["n"]/df["M"])

#depleted
df["odds_ratio"] = (df["n"]/df["M"]) / (df["k"]/df["N"])
