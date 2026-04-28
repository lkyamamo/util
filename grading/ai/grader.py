#!/usr/bin/env python3

import os
import re
import json
import argparse
import pdfplumber
import pandas as pd
import nbformat
from tqdm import tqdm
from openai import OpenAI

MODEL = "gpt-5.4"
client = OpenAI()


RUBRIC_GENERATION_PROMPT = """
You are converting a course assignment sheet into a grading rubric for an automated compliance checker.

Your job is to translate the sheet into STRICT JSON with this exact schema:
{
  "requirements": [
    {
      "id": "P1_a_code_shape",
      "type": "code",
      "description": "Problem 1(a): Reports the shape of X and y."
    }
  ]
}

Allowed requirement types:
- "code": verify from code only
- "written": verify from written text only
- "integration": verify that code-produced results are also discussed or presented in the written submission

Rules:
- Include every explicit student-facing task that produces a checkable deliverable.
- Split multipart tasks into separate requirements when they can be checked independently.
- Preserve problem/part structure in the IDs and descriptions.
- Use short stable IDs like P1_a_code_shape, P1_d_written_f1_interpretation, P2_e_integration_plot.
- Descriptions must be concrete, checkable, and phrased as submission requirements.
- Prefer one requirement per distinct deliverable, computation, plot, table, or explanation.
- For plots/figures/tables required by the sheet, create:
  - a "code" requirement for generating it
  - an "integration" requirement only when the sheet expects it to appear in the submitted report or be discussed there
- Use "written" for explicit explanation / interpretation / discussion prompts.
- Do not include course logistics such as due dates, compression format, upload platform, or other administrative instructions unless they are part of the academic deliverables themselves.
- Do not include background theory unless students are explicitly asked to explain or interpret it.
- Do not invent requirements beyond what the sheet asks.
- Return STRICT JSON only. No markdown. No commentary.
"""


# ----------------------------
# FILE EXTRACTION
# ----------------------------

def extract_pdf_text(pdf_path):
    text = ""
    try:
        with pdfplumber.open(pdf_path) as pdf:
            for page in pdf.pages:
                text += page.extract_text() or ""
    except Exception as e:
        print(f"PDF extraction failed: {pdf_path} ({e})")
    return text


def extract_ipynb_code(ipynb_path):
    code_text = ""
    try:
        nb = nbformat.read(ipynb_path, as_version=4)
        for cell in nb.cells:
            if cell.cell_type == "code":
                code_text += cell.source + "\n\n"
    except Exception as e:
        print(f"Notebook extraction failed: {ipynb_path} ({e})")
    return code_text


def collect_submission_content(student_path):
    written_text = ""
    code_text = ""

    for root, _, files in os.walk(student_path):
        for file in files:
            full_path = os.path.join(root, file)

            if file.endswith(".pdf"):
                written_text += extract_pdf_text(full_path) + "\n"

            elif file.endswith(".py"):
                with open(full_path, "r", encoding="utf-8", errors="ignore") as f:
                    code_text += f.read() + "\n\n"

            elif file.endswith(".ipynb"):
                code_text += extract_ipynb_code(full_path)

    return written_text.strip(), code_text.strip()


# ----------------------------
# AI EVALUATION
# ----------------------------

def normalize_requirements(rubric):
    """
    Support two rubric shapes:
    1) {"requirements": [{"id","type","description"}, ...]}
    2) Flat lists: {"written_requirements":[...], "code_requirements":[...], "integration_requirements":[...]}
       In this case, we use the description string as the id to preserve compatibility with downstream tooling.
    """
    if isinstance(rubric, dict) and isinstance(rubric.get("requirements"), list):
        return rubric["requirements"]

    reqs = []
    key_to_type = {
        "written_requirements": "written",
        "code_requirements": "code",
        "integration_requirements": "integration",
    }
    if not isinstance(rubric, dict):
        return reqs

    for key, rtype in key_to_type.items():
        items = rubric.get(key, [])
        if not isinstance(items, list):
            continue
        for desc in items:
            if not isinstance(desc, str):
                continue
            reqs.append({"id": desc, "type": rtype, "description": desc})
    return reqs


def _extract_json_from_response(text):
    """Parse JSON from model response, tolerating markdown code fences and stray newlines."""
    s = (text or "").strip()
    # Strip markdown code block if present
    if s.startswith("```"):
        lines = s.split("\n")
        if lines[0].strip().startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].strip() == "```":
            lines = lines[:-1]
        s = "\n".join(lines)

    def _try_parse(raw):
        # Try as-is first
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            pass
        # Remove literal newlines/tabs that the model may have injected mid-token
        cleaned = re.sub(r'[\n\r\t]+', ' ', raw)
        try:
            return json.loads(cleaned)
        except json.JSONDecodeError:
            return None

    result = _try_parse(s)
    if result is not None:
        return result

    # Fallback: extract the first JSON object substring
    start = s.find("{")
    end = s.rfind("}")
    if start != -1 and end != -1 and end > start:
        return _try_parse(s[start : end + 1])
    return None


def _normalize_rubric_keys(rubric):
    """Strip whitespace from all keys and string values in each requirement dict."""
    reqs = rubric.get("requirements")
    if not isinstance(reqs, list):
        return rubric
    cleaned = []
    for req in reqs:
        if not isinstance(req, dict):
            cleaned.append(req)
            continue
        cleaned.append({k.strip(): (v.strip() if isinstance(v, str) else v)
                        for k, v in req.items()})
    return {**rubric, "requirements": cleaned}


def _validate_generated_rubric(rubric):
    rubric = _normalize_rubric_keys(rubric)
    requirements = normalize_requirements(rubric)
    if not requirements:
        raise ValueError(
            "Generated rubric has no requirements. "
            "Expected 'requirements' or compatible flat requirement lists."
        )

    valid_types = {"code", "written", "integration"}
    for req in requirements:
        if not isinstance(req, dict):
            raise ValueError("Generated rubric contains a non-dict requirement.")
        if not req.get("id") or not req.get("description"):
            raise ValueError("Each generated requirement must have 'id' and 'description'.")
        if req.get("type") not in valid_types:
            raise ValueError(
                f"Invalid generated requirement type: {req.get('type')!r}. "
                "Expected one of: code, written, integration."
            )
    return requirements


def _chat(messages):
    """Call the model, retrying without temperature if the model doesn't support it."""
    try:
        response = client.chat.completions.create(
            model=MODEL,
            temperature=0,
            messages=messages,
        )
    except Exception as e:
        if "temperature" in str(e).lower() or "unsupported" in str(e).lower():
            response = client.chat.completions.create(
                model=MODEL,
                messages=messages,
            )
        else:
            raise
    return response.choices[0].message.content


def generate_rubric_from_assignment_text(assignment_text):
    prompt = f"""
{RUBRIC_GENERATION_PROMPT}

Assignment sheet text:
\"\"\"
{assignment_text[:30000]}
\"\"\"
"""

    content = _chat([{"role": "user", "content": prompt}])
    rubric = _extract_json_from_response(content)
    if rubric is None:
        print(f"\n--- Raw model response (rubric generation) ---\n{content}\n---\n")
        raise ValueError("Could not parse rubric JSON from assignment-sheet generation response.")
    rubric = _normalize_rubric_keys(rubric)
    _validate_generated_rubric(rubric)
    return rubric


def load_or_generate_rubric(args):
    if args.rubric:
        with open(args.rubric, "r") as f:
            rubric = json.load(f)
        requirements = normalize_requirements(rubric)
        if not requirements:
            raise ValueError(
                "Rubric has no requirements. Expected either 'requirements' or "
                "'written_requirements'/'code_requirements'/'integration_requirements'."
            )
        return rubric, requirements

    if args.rubric_output and os.path.exists(args.rubric_output) and not args.regrade_all:
        print(f"Using existing rubric: {args.rubric_output}")
        with open(args.rubric_output, "r") as f:
            rubric = json.load(f)
        requirements = normalize_requirements(rubric)
        if requirements:
            return rubric, requirements
        print("Existing rubric has no valid requirements — regenerating.")

    if not args.assignment_sheet:
        raise ValueError("Provide either --rubric or --assignment_sheet.")

    assignment_text = extract_pdf_text(args.assignment_sheet)
    if not assignment_text.strip():
        raise ValueError(f"No text could be extracted from assignment sheet: {args.assignment_sheet}")

    rubric = generate_rubric_from_assignment_text(assignment_text)
    requirements = normalize_requirements(rubric)

    if args.rubric_output:
        with open(args.rubric_output, "w", encoding="utf-8") as f:
            json.dump(rubric, f, indent=2)
        print(f"Saved generated rubric to {args.rubric_output}")

    return rubric, requirements


def evaluate_requirements(requirements, written_text, code_text):

    prompt = f"""
You are an assignment compliance evaluator.

Rules:
- For type "code": evaluate using code only.
- For type "written": evaluate using written text only.
- For type "integration": evaluate consistency between written and code.
- Be lenient for written.
- Mark NOT_MET only when the requirement is clearly missing or violated.
- Prefer MET when the code or written text clearly satisfies the requirement. Use UNCERTAIN only when you truly cannot verify (e.g. content is missing, unreadable, or genuinely ambiguous). Do not use UNCERTAIN when the requirement is satisfied.

You MUST return one result per requirement. For each result, copy the exact "id" string from the requirement—do not modify or abbreviate it.

Return STRICT JSON only (no markdown):

{{
  "results": [
    {{
      "id": "<exact id from requirements>",
      "status": "MET | NOT_MET | UNCERTAIN"
    }}
  ]
}}

Requirements:
{json.dumps(requirements, indent=2)}

Written:
\"\"\"
{written_text[:12000]}
\"\"\"

Code:
\"\"\"
{code_text[:25000]}
\"\"\"
"""

    return _chat([{"role": "user", "content": prompt}])


# ----------------------------
# MAIN
# ----------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--submissions")
    parser.add_argument("--rubric")
    parser.add_argument("--assignment_sheet")
    parser.add_argument("--rubric_output")
    parser.add_argument("--output_prefix", default="results")
    parser.add_argument("--regrade-all", action="store_true",
                        help="Re-grade all students even if results already exist")

    args = parser.parse_args()
    rubric, requirements = load_or_generate_rubric(args)

    if not args.submissions:
        print("Rubric generation complete.")
        return

    code_requirements = [r for r in requirements if r["type"] == "code"]
    written_requirements = [
        r for r in requirements if r["type"] in ["written", "integration"]
    ]

    students = [
        s for s in os.listdir(args.submissions)
        if os.path.isdir(os.path.join(args.submissions, s))
    ]

    existing_students = set()
    code_rows = []
    written_rows = []
    if not args.regrade_all:
        code_csv_path = f"{args.output_prefix}_code.csv"
        written_csv_path = f"{args.output_prefix}_written.csv"
        if os.path.exists(code_csv_path):
            existing_df = pd.read_csv(code_csv_path)
            existing_students = set(existing_df["Student"].astype(str).tolist())
            code_rows = existing_df.to_dict("records")
        if os.path.exists(written_csv_path):
            written_rows = pd.read_csv(written_csv_path).to_dict("records")
        if existing_students:
            print(f"Skipping {len(existing_students)} already-graded student(s).")

    students_to_grade = [s for s in students if s not in existing_students]
    print(f"Grading {len(students_to_grade)} student(s).")

    for student in tqdm(students_to_grade):

        student_path = os.path.join(args.submissions, student)
        written_text, code_text = collect_submission_content(student_path)

        # Only evaluate code requirements with AI; written/integration are flagged for manual review.
        ai_response = evaluate_requirements(
            code_requirements, written_text, code_text
        )

        parsed = []
        obj = _extract_json_from_response(ai_response)
        if isinstance(obj, dict) and isinstance(obj.get("results"), list):
            parsed = obj["results"]

        # Build result map by id; fall back to index when lengths match
        result_map = {}
        for r in parsed:
            if isinstance(r, dict) and "id" in r and "status" in r:
                result_map[r["id"]] = r["status"]
        if len(parsed) == len(code_requirements):
            for i, req in enumerate(code_requirements):
                req_id = req["id"]
                if req_id not in result_map and isinstance(parsed[i], dict):
                    result_map[req_id] = parsed[i].get("status", "UNCERTAIN")

        def status_for(req):
            return result_map.get(req["id"], "UNCERTAIN")

        # Code CSV
        code_row = {"Student": student}
        for req in code_requirements:
            code_row[req["id"]] = status_for(req)
        code_rows.append(code_row)

        # Written CSV — flagged for manual review, not AI-graded
        written_row = {"Student": student}
        for req in written_requirements:
            written_row[req["id"]] = "NEEDS_REVIEW"
        written_rows.append(written_row)

    pd.DataFrame(code_rows).to_csv(
        f"{args.output_prefix}_code.csv", index=False
    )
    pd.DataFrame(written_rows).to_csv(
        f"{args.output_prefix}_written.csv", index=False
    )

    print("Grading complete.")
    print(f"Generated:")
    print(f"- {args.output_prefix}_code.csv")
    print(f"- {args.output_prefix}_written.csv")


if __name__ == "__main__":
    main()