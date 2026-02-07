# Agent: Data Downloader

## Purpose
Downloads biological databases from public sources and prepares them for adapter processing.

## Guidelines
- Download to template_package/data/<database_name>/
- Create download manifest JSON with file checksums
- Create analysis report in ANALYSIS.md
- Use streaming downloads for large files
- Support resume on interrupted downloads
- Respect rate limits on APIs
