#!/usr/bin/env python3
"""
Script to create a new release by updating the version in pyproject.toml.
This script automates semantic version bumping and safely commits/pushes only
pyproject.toml.
"""

import argparse
import re
import subprocess
import sys

BUMP_LEVELS = {"none": 0, "patch": 1, "minor": 2, "major": 3}
LEVEL_TO_BUMP = {0: "none", 1: "patch", 2: "minor", 3: "major"}
PYPROJECT_FILE = "pyproject.toml"


def get_current_version():
    """Get the current version from pyproject.toml"""
    try:
        with open(PYPROJECT_FILE, "r", encoding="utf-8") as f:
            content = f.read()
        version = _extract_project_version(content)
        if not version:
            raise ValueError("Could not find [project].version in pyproject.toml")
        return version
    except (OSError, ValueError) as e:
        print(f"Error reading pyproject.toml: {e}")
        sys.exit(1)


def update_version(new_version):
    """Update only the [project] version line in pyproject.toml."""
    try:
        with open(PYPROJECT_FILE, "r", encoding="utf-8") as f:
            content = f.read()

        updated, changed = _replace_project_version(content, new_version)
        if not changed:
            print("No version change applied to pyproject.toml")
            return False

        with open(PYPROJECT_FILE, "w", encoding="utf-8", newline="") as f:
            f.write(updated)

        print(f"Updated version to {new_version} in pyproject.toml")
        return True
    except (OSError, ValueError) as e:
        print(f"Error updating pyproject.toml: {e}")
        return False


def _extract_project_version(content):
    """Extract [project] version from pyproject content."""
    in_project = False
    for line in content.splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_project = stripped == "[project]"
            continue
        if in_project and re.match(r'^version\s*=\s*"([^"]+)"\s*$', stripped):
            return re.match(r'^version\s*=\s*"([^"]+)"\s*$', stripped).group(1)
    return None


def _replace_project_version(content, new_version):
    """Replace only [project] version line, preserving formatting elsewhere."""
    lines = content.splitlines(keepends=True)
    in_project = False
    changed = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_project = stripped == "[project]"
            continue
        if in_project:
            m = re.match(r'^(\s*version\s*=\s*")([^"]+)(".*?)(\r?\n?)$', line)
            if m:
                old_version = m.group(2)
                if old_version != new_version:
                    lines[i] = f'{m.group(1)}{new_version}{m.group(3)}{m.group(4)}'
                    changed = True
                return "".join(lines), changed
    raise ValueError("Could not locate [project] version line in pyproject.toml")


def validate_version(version):
    """Validate that the version follows semantic versioning"""
    pattern = r'^\d+\.\d+\.\d+$'
    if not re.match(pattern, version):
        print(f"Error: Version '{version}' must follow semantic versioning (e.g., 1.0.0)")
        return False
    return True


def parse_version(version):
    """Parse semantic version string into a tuple of ints."""
    if not validate_version(version):
        raise ValueError(f"Invalid semantic version: {version}")
    major, minor, patch = version.split(".")
    return int(major), int(minor), int(patch)


def bump_version(current_version, bump):
    """Compute next semantic version from current version and bump type."""
    major, minor, patch = parse_version(current_version)
    if bump == "major":
        return f"{major + 1}.0.0"
    if bump == "minor":
        return f"{major}.{minor + 1}.0"
    if bump == "patch":
        return f"{major}.{minor}.{patch + 1}"
    raise ValueError(f"Unsupported bump type: {bump}")


def run_command(cmd, check=True):
    """Run a shell command"""
    try:
        result = subprocess.run(cmd, shell=True, check=check, capture_output=True, text=True)
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        return False, e.stdout, e.stderr


def has_pyproject_version_change():
    """Return True if pyproject.toml has unstaged/staged changes."""
    success, stdout, stderr = run_command("git diff -- pyproject.toml", check=False)
    if not success:
        raise RuntimeError(f"Unable to check pyproject.toml diff: {stderr.strip()}")
    success_cached, stdout_cached, stderr_cached = run_command(
        "git diff --cached -- pyproject.toml", check=False
    )
    if not success_cached:
        raise RuntimeError(f"Unable to check staged pyproject.toml diff: {stderr_cached.strip()}")
    return bool(stdout.strip() or stdout_cached.strip())


def get_last_tag():
    """Return the most recent git tag or None if repository has no tags."""
    success, stdout, _ = run_command("git describe --tags --abbrev=0", check=False)
    if not success:
        return None
    tag = stdout.strip()
    return tag or None


def get_current_branch():
    """Return current git branch name."""
    success, stdout, stderr = run_command("git rev-parse --abbrev-ref HEAD", check=False)
    if not success:
        raise RuntimeError(f"Unable to detect current branch: {stderr.strip()}")
    branch = stdout.strip()
    if not branch:
        raise RuntimeError("Unable to detect current branch: empty branch name")
    return branch


def get_commits_and_files_since_ref(git_ref):
    """Collect commit messages and changed files since git_ref."""
    if git_ref:
        commit_range = f"{git_ref}..HEAD"
        diff_range = f"{git_ref}..HEAD"
    else:
        commit_range = "--all"
        diff_range = "HEAD"

    success, commits_out, commits_err = run_command(
        f'git log {commit_range} --pretty=format:"%H%x1f%s%x1f%b%x1e"',
        check=False
    )
    if not success:
        raise RuntimeError(f"Unable to read git log: {commits_err.strip()}")

    success, files_out, files_err = run_command(
        f"git diff --name-only {diff_range}",
        check=False
    )
    if not success:
        raise RuntimeError(f"Unable to read git diff: {files_err.strip()}")

    commits = []
    for entry in commits_out.strip("\n\x1e").split("\x1e"):
        if not entry.strip():
            continue
        parts = entry.split("\x1f")
        if len(parts) < 3:
            continue
        commits.append({
            "hash": parts[0].strip(),
            "subject": parts[1].strip(),
            "body": parts[2].strip(),
        })

    files = [line.strip() for line in files_out.splitlines() if line.strip()]
    return commits, files


def _classify_files(changed_files):
    """Classify changed files for reporting and fallback bump logic."""
    categories = {
        "source": [],
        "tests": [],
        "docs": [],
        "scripts": [],
        "config": [],
        "other": [],
    }
    for path in changed_files:
        norm = path.replace("\\", "/")
        if norm.startswith("src/pycequeau/"):
            categories["source"].append(path)
        elif norm.startswith("test") or norm.startswith("tests/"):
            categories["tests"].append(path)
        elif norm.startswith("docs/") or norm.lower().startswith("readme"):
            categories["docs"].append(path)
        elif norm.startswith("scripts/"):
            categories["scripts"].append(path)
        elif norm in {"pyproject.toml", "MANIFEST.in", "environment.yml", "environment_pyceq.yml"}:
            categories["config"].append(path)
        else:
            categories["other"].append(path)
    return categories


def detect_bump_type():
    """Detect the recommended semantic version bump from commit metadata and changed files."""
    last_tag = get_last_tag()
    commits, files = get_commits_and_files_since_ref(last_tag)
    categories = _classify_files(files)

    level = BUMP_LEVELS["none"]
    reasons = []

    breaking_re = re.compile(r"^[a-zA-Z]+(\([^)]+\))?!:")
    feat_re = re.compile(r"^feat(\([^)]+\))?:")
    fix_re = re.compile(r"^(fix|perf|refactor)(\([^)]+\))?:")

    for commit in commits:
        subject = commit["subject"]
        body = commit["body"]
        short = commit["hash"][:7]
        if "BREAKING CHANGE" in body or breaking_re.match(subject):
            level = max(level, BUMP_LEVELS["major"])
            reasons.append(f"{short}: breaking change detected in commit message")
        elif feat_re.match(subject):
            level = max(level, BUMP_LEVELS["minor"])
            reasons.append(f"{short}: feature commit detected ({subject})")
        elif fix_re.match(subject):
            level = max(level, BUMP_LEVELS["patch"])
            reasons.append(f"{short}: fix/perf/refactor commit detected ({subject})")

    # Fallback heuristics when commit messages are not conventional
    if level == BUMP_LEVELS["none"]:
        if categories["source"]:
            level = BUMP_LEVELS["patch"]
            reasons.append("source code changes detected (fallback to patch)")
        elif categories["scripts"] or categories["config"]:
            level = BUMP_LEVELS["patch"]
            reasons.append("release/build/config changes detected (fallback to patch)")
        elif categories["docs"] and not (categories["tests"] or categories["other"]):
            reasons.append("docs-only changes detected")
        elif categories["tests"]:
            level = BUMP_LEVELS["patch"]
            reasons.append("test-only changes detected (fallback to patch)")
        elif categories["other"]:
            level = BUMP_LEVELS["patch"]
            reasons.append("unclassified changes detected (fallback to patch)")

    return {
        "bump": LEVEL_TO_BUMP[level],
        "last_tag": last_tag,
        "reasons": reasons,
        "categories": categories,
        "commit_count": len(commits),
    }


def print_detection_summary(detection):
    """Print auto-detection result and rationale."""
    print("\nAuto-detected semantic version bump")
    print("----------------------------------")
    print(f"Last tag: {detection['last_tag'] or 'None'}")
    print(f"Commits inspected: {detection['commit_count']}")
    print(f"Detected bump: {detection['bump']}")
    print("Changed files by category:")
    for category, files in detection["categories"].items():
        print(f"  - {category}: {len(files)}")
    if detection["reasons"]:
        print("Reasons:")
        for reason in detection["reasons"]:
            print(f"  - {reason}")
    else:
        print("Reasons:")
        print("  - No changes requiring a version bump were detected")


def confirm(prompt):
    """Prompt for yes/no confirmation."""
    return input(prompt).strip().lower() == "y"


def resolve_new_version(args, current_version):
    """Resolve target version from manual input or bump strategy."""
    if args.version and args.bump:
        print("Error: Provide either a manual version OR --bump, not both.")
        sys.exit(1)

    if not args.version and not args.bump:
        args.bump = "auto"

    if args.version:
        if not validate_version(args.version):
            sys.exit(1)
        print(f"New version: {args.version} (manual)")
        return args.version

    bump = args.bump
    if args.bump == "auto":
        try:
            detection = detect_bump_type()
        except RuntimeError as exc:
            print(f"Error during auto-detection: {exc}")
            sys.exit(1)

        print_detection_summary(detection)
        if detection["bump"] == "none":
            print("\nNo version bump required from detected changes.")
            if not confirm("Proceed with patch bump anyway? (y/N): "):
                print("Aborted")
                sys.exit(0)
            bump = "patch"
        else:
            bump = detection["bump"]

    new_version = bump_version(current_version, bump)
    print(f"New version: {new_version} (--bump {args.bump})")
    return new_version


def resolve_target_branch(push_branch):
    """Resolve branch to push to."""
    try:
        current_branch = get_current_branch()
    except RuntimeError as exc:
        print(f"Error: {exc}")
        sys.exit(1)

    target_branch = current_branch if push_branch == "auto" else push_branch
    if target_branch == "main":
        print("Warning: target push branch is 'main'. If protected, push will fail.")
    return current_branch, target_branch


def commit_and_push(new_version, target_branch):
    """Commit and push only pyproject.toml."""
    if not confirm(f"\nConfirm commit and push to origin/{target_branch}? (y/N): "):
        print("Aborted before commit/push")
        sys.exit(0)

    print("Committing version change...")
    success, _, stderr = run_command("git add pyproject.toml")
    if not success:
        print(f"Error adding pyproject.toml: {stderr}")
        sys.exit(1)

    success, _, stderr = run_command(
        f'git commit -m "Bump version to {new_version}" -- pyproject.toml'
    )
    if not success:
        print(f"Error committing changes: {stderr}")
        sys.exit(1)
    print("Changes committed successfully")

    print(f"Pushing commit to origin/{target_branch} ...")
    success, _, stderr = run_command(f"git push origin {target_branch}")
    if not success:
        print(f"Error pushing changes: {stderr}")
        sys.exit(1)

    print("✅ Version bump commit created and pushed for PR workflow.")
    print(
        "After PR merge to main, GitHub Actions will create tag/release "
        "and publish from pyproject version."
    )


def main():
    """Main function to handle command line arguments and execute release process."""
    parser = argparse.ArgumentParser(description="Create a new release")
    parser.add_argument(
        "version",
        nargs="?",
        help="Manual new version number (e.g., 1.0.0). If omitted, use --bump.",
    )
    parser.add_argument(
        "--bump",
        choices=["major", "minor", "patch", "auto"],
        help="Compute the next version from current version.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes",
    )
    parser.add_argument(
        "--skip-commit",
        action="store_true",
        help="Skip committing changes (useful for testing)",
    )
    parser.add_argument(
        "--push-branch",
        default="auto",
        help="Branch to push. Default: current branch (auto)."
    )

    args = parser.parse_args()

    current_version = get_current_version()
    print(f"Current version: {current_version}")
    new_version = resolve_new_version(args, current_version)

    if current_version == new_version:
        print("Warning: Version is already set to the requested version")
        if not confirm("Continue anyway? (y/N): "):
            print("Aborted")
            sys.exit(0)

    current_branch, target_branch = resolve_target_branch(args.push_branch)

    if args.dry_run:
        print("\nDry run - would perform the following actions:")
        print(f"1. Update pyproject.toml version from {current_version} to {new_version}")
        print(f"2. Commit changes with message: 'Bump version to {new_version}'")
        print(f"3. Push commit to origin/{target_branch} (current branch: {current_branch})")
        return

    # Update version in pyproject.toml
    if not update_version(new_version):
        sys.exit(1)

    try:
        if not has_pyproject_version_change():
            print("Error: pyproject.toml version was not changed. Aborting.")
            sys.exit(1)
    except RuntimeError as exc:
        print(f"Error: {exc}")
        sys.exit(1)

    if args.skip_commit:
        print("Skipped commit and tag creation (--skip-commit flag used)")
        print("You can manually commit and push when ready")
        return

    commit_and_push(new_version, target_branch)


if __name__ == "__main__":
    main()
