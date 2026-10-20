from pathlib import Path
import pandas as pd
TRACE_FILENAMES = [
    "https://github.com/Azure/AzurePublicDataset/releases/download/dataset-llm-2024/AzureLLMInferenceTrace_code_1week.csv",
]
tokens = pd.read_csv(TRACE_FILENAMES[0], parse_dates=["TIMESTAMP"], date_format="mixed", dayfirst=False)
tokens["bin"] = pd.to_datetime(tokens.pop("TIMESTAMP"), format="%Y-%m-%d %H:%M:%S.%f%z",
    utc=True).dt.floor("5min")
means = (tokens.groupby("bin", as_index=False)[["GeneratedTokens"]]
    .mean()
    .rename(columns={"GeneratedTokens": "m_gen"})
)
means.to_csv("training_loss/InferenceTrace_code_5min.csv", index=False)
