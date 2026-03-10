import json
import re
from pathlib import Path
from collections import defaultdict

KEYWORDS = {
    "P1_R1": ["MSE", "MAE"],
    "P1_R2": ["scatter", "y = x"],
    "P1_R3": ["outlier", "robust"],
    "P1_R4": ["RMSE"],
    "P2_R4": ["15", "23", "30"],
    "P3_R1": ["np.random.seed(42)", "200"],
    "P3_R5": ["correlation", "Pearson"]
}

FIGURE_EXTENSIONS = {".png", ".jpg", ".jpeg", ".pdf"}

def extract_evidence(submission):
    evidence = defaultdict(list)

    contents = submission["contents"]
    files = submission["file_tree"]

    # Detect figures
    for f in files:
        if Path(f).suffix.lower() in FIGURE_EXTENSIONS:
            evidence["figures"].append(f)

    # Search text for keywords
    for req_id, words in KEYWORDS.items():
        for file, text in contents.items():
            for w in words:
                if w.lower() in text.lower():
                    snippet = extract_snippet(text, w)
                    evidence[req_id].append(f"{file}: ...{snippet}...")

    return dict(evidence)

def extract_snippet(text, keyword, window=40):
    idx = text.lower().find(keyword.lower())
    if idx == -1:
        return ""
    start = max(0, idx - window)
    end = min(len(text), idx + window)
    return text[start:end]
