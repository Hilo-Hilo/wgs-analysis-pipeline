# GitHub Upload Instructions - PRIVATE Repository First

## ⚠️ IMPORTANT: Create as PRIVATE Repository First

Follow these steps to safely upload your WGS pipeline to GitHub as a **PRIVATE** repository first. You can make it public later after final review.

## Step 1: Create PRIVATE GitHub Repository

### Via GitHub Website:
1. Go to https://github.com/new
2. Repository name: `wgs-analysis-pipeline`
3. Description: "Whole genome sequencing analysis pipeline and tools"
4. **IMPORTANT: Select "Private" (not Public)**
5. DO NOT initialize with README (we already have one)
6. DO NOT add .gitignore (we already have one)
7. DO NOT choose a license (we already have MIT license)
8. Click "Create repository"

## Step 2: Review Files Locally

Before pushing, double-check that NO personal data is included:

```bash
cd ~/wgs-analysis-public

# Check for any personal identifiers
grep -r "E200032534" . 2>/dev/null
grep -r "hansonwen" . 2>/dev/null

# Check file sizes (no large genomic files should be present)
find . -type f -size +10M -exec ls -lh {} \;

# Review the files that will be committed
git status

# Check what's being ignored
git status --ignored
```

## Step 3: Commit and Push to PRIVATE Repository

```bash
# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: WGS analysis pipeline without personal data"

# Add your private GitHub repository as remote
# Replace YOUR_USERNAME with your GitHub username
git remote add origin https://github.com/YOUR_USERNAME/wgs-analysis-pipeline.git

# Push to main branch
git branch -M main
git push -u origin main
```

## Step 4: Verify on GitHub (While Still Private)

1. Go to your repository: https://github.com/YOUR_USERNAME/wgs-analysis-pipeline
2. Verify it shows "Private" label
3. Check all directories and files
4. Ensure NO genomic data files (.bam, .vcf, .fastq) are present
5. Verify scripts have SAMPLE001 instead of personal IDs

## Step 5: Final Safety Checks

### Check for Sensitive Information:
- [ ] No .vcf, .vcf.gz files
- [ ] No .bam, .sam files  
- [ ] No .fastq, .fq files
- [ ] No files with "E200032534" in name
- [ ] No files with personal names
- [ ] No medical reports with real data
- [ ] No variant lists with actual mutations
- [ ] Scripts use "SAMPLE001" not real IDs
- [ ] Documentation is anonymized

### Repository Contains Only:
- [x] Shell scripts (.sh)
- [x] Python scripts (.py)
- [x] Documentation (.md)
- [x] Templates (no real data)
- [x] Configuration examples
- [x] .gitignore
- [x] LICENSE
- [x] README

## Step 6: Optional - Add Collaborators for Review

While private, you can add trusted collaborators to review:

1. Go to Settings → Manage access
2. Click "Invite a collaborator"
3. Add their GitHub username
4. They can review and confirm no personal data

## Step 7: Make Public (When Ready)

**Only after thorough review:**

1. Go to Settings (bottom of sidebar)
2. Scroll to "Danger Zone"
3. Click "Change visibility"
4. Select "Make public"
5. Type repository name to confirm
6. Click "I understand, make this repository public"

## Alternative: Use GitHub CLI

If you have GitHub CLI installed:

```bash
# Create private repo
gh repo create wgs-analysis-pipeline --private --source=. --remote=origin

# Push code
git push -u origin main

# Later, make public (after review)
gh repo edit YOUR_USERNAME/wgs-analysis-pipeline --visibility public
```

## What's Been Anonymized

- Sample ID: E200032534 → SAMPLE001
- User paths: /Users/hansonwen → /Users/user
- Personal reports: Removed entirely
- Variant data: Excluded via .gitignore
- Medical findings: Converted to templates

## Additional Safety Tips

1. **Use GitHub's Secret Scanning**: GitHub automatically scans for exposed secrets
2. **Enable Branch Protection**: Prevent accidental pushes to main
3. **Review GitHub Insights**: Check what files are being accessed
4. **Use `.gitattributes`**: Mark certain files as binary to prevent diffing

## If You Accidentally Push Personal Data

If you accidentally commit personal data:

```bash
# DO NOT PANIC - If repo is private, you have time

# Remove file from history (example)
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch PATH_TO_FILE" \
  --prune-empty --tag-name-filter cat -- --all

# Force push
git push origin --force --all

# Or use BFG Repo-Cleaner (easier)
# Download from: https://rtyley.github.io/bfg-repo-cleaner/
java -jar bfg.jar --delete-files FILENAME
git push --force
```

## Final Checklist Before Making Public

- [ ] Reviewed all files in GitHub web interface
- [ ] Checked commit history for sensitive data
- [ ] Confirmed .gitignore is working properly
- [ ] Tested that others can't see private data
- [ ] Documentation is helpful but anonymous
- [ ] Ready to share with the community

## Support

If you're unsure about anything, keep the repository private and seek advice from:
- GitHub Support
- Bioinformatics communities
- Privacy-focused developers

Remember: **It's always better to be overly cautious with genomic data!**