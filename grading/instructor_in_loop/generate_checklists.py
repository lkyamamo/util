import json
import yaml
from pathlib import Path
from evidence_extractor import extract_evidence

SUBMISSIONS_DIR = "prepared_submissions"
REQUIREMENTS_FILE = "requirements.yaml"
OUTPUT_DIR = "checklists"

Path(OUTPUT_DIR).mkdir(exist_ok=True)

def main():
    requirements = yaml.safe_load(open(REQUIREMENTS_FILE))

    for submission_file in Path(SUBMISSIONS_DIR).glob("*.json"):
        submission = json.load(open(submission_file))
        student = submission["student"]

        evidence = extract_evidence(submission)

        lines = []
        lines.append(f"# PHYS 515 – Exercise Sheet 2")
        lines.append(f"## Student: {student}\n")

        lines.append("### File Tree")
        for f in submission["file_tree"]:
            lines.append(f"- {f}")

        lines.append("\n---\n")

        for section, reqs in requirements.items():
            if section == "assignment_metadata":
                continue

            lines.append(f"### {section.replace('_', ' ').title()}\n")

            for r in reqs:
                rid = r["id"]
                lines.append(f"[ ] {rid} {r['description']}")
                lines.append("    Evidence found:")

                if rid in evidence:
                    for e in evidence[rid]:
                        lines.append(f"    - {e}")
                else:
                    lines.append("    - (none detected)")

                lines.append("")

        out = Path(OUTPUT_DIR) / f"{student}.md"
        out.write_text("\n".join(lines))

if __name__ == "__main__":
    main()
