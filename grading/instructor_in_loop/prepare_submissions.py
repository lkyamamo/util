import os
import json
import subprocess
from pathlib import Path

TEXT_EXTENSIONS = {".txt", ".md", ".py", ".ipynb"}
PDF_EXTENSIONS = {".pdf"}

def extract_pdf_text(pdf_path):
    try:
        result = subprocess.run(
            ["pdftotext", str(pdf_path), "-"],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except Exception as e:
        return f"[PDF extraction failed: {e}]"

def build_file_tree(root):
    tree = []
    for path in sorted(root.rglob("*")):
        if path.is_file():
            tree.append(str(path.relative_to(root)))
    return tree

def extract_contents(root):
    contents = {}
    for path in root.rglob("*"):
        if path.suffix.lower() in TEXT_EXTENSIONS:
            try:
                contents[str(path.relative_to(root))] = path.read_text(errors="ignore")
            except Exception as e:
                contents[str(path.relative_to(root))] = f"[Read error: {e}]"
        elif path.suffix.lower() in PDF_EXTENSIONS:
            contents[str(path.relative_to(root))] = extract_pdf_text(path)
    return contents

def process_submission(student_dir, output_dir):
    student_dir = Path(student_dir)
    output = {
        "student": student_dir.name,
        "file_tree": build_file_tree(student_dir),
        "contents": extract_contents(student_dir)
    }

    output_path = Path(output_dir) / f"{student_dir.name}.json"
    output_path.write_text(json.dumps(output, indent=2))

def main(submissions_root, output_dir):
    Path(output_dir).mkdir(exist_ok=True)
    for student_dir in Path(submissions_root).iterdir():
        if student_dir.is_dir():
            process_submission(student_dir, output_dir)

if __name__ == "__main__":
    main("submissions/", "prepared_submissions/")
