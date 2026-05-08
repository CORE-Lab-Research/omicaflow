# Contributing to OmicaFlow

Thank you for your interest in contributing to OmicaFlow! This guide outlines how to contribute to the project.

## Code of Conduct
By participating, you agree to uphold our lab's collaborative and respectful environment.

## How to Contribute
1. **Fork the repository** to your own GitHub account
2. **Create a feature branch** (`git checkout -b feature/your-feature-name`)
3. **Make your changes** following the coding standards below
4. **Commit clearly** with descriptive messages (see Commit Guidelines)
5. **Push to your fork** and open a Pull Request to the `main` branch

## Coding Standards
- **R Scripts**: Follow Tidyverse style guide; no unnecessary comments
- **Python Scripts**: Follow PEP8; use type hints where possible
- **Snakemake Rules**: Use clear input/output wildcards; pin tool versions
- **Config**: Only modify `config/base.yaml` for adjustable parameters

## Commit Guidelines
- Use imperative mood ("Add module" not "Added module")
- Keep commits atomic (one logical change per commit)
- Reference issue/task numbers if applicable (e.g., "T006: Implement QC module")

## Development Workflow
- New modules require:
  1. New directory under `modules/`
  2. Corresponding Snakemake rule in `workflow/rules/`
  3. Toggle flag in `config/base.yaml`
  4. Update `OMICAFLOW_MVP_PLAN.md` progress table

## Testing
- Run `snakemake -n` for dry-run validation before committing
- Ensure new modules work with the existing pipeline structure

## Questions?
Open an issue on the [GitHub repo](https://github.com/CORE-Lab-Research/omicaflow/issues) for discussion.