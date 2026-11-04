# NERDSS user guide (online edition)

This folder holds the NERDSS user guide as Markdown, plus the web page built from it.

| File | What it is |
|------|------------|
| `nerdss_user_guide.md` | **The source. Edit this file.** Readable as is on GitHub. |
| `nerdss_user_guide.html` | The page on the website. Generated, so do not edit it by hand. |
| `media/` | Figures used by the guide (taken from the Word document). |
| `build.sh` | Regenerates the table of contents in the `.md` and rebuilds the `.html`. |
| `template.html` | Pandoc template: the website's header, menu, footer and scripts around the guide. |
| `filters.lua` | Pandoc filter: title from the first heading, scrollable tables, warning callouts. |
| `toc.tpl` | Tiny pandoc template used by `build.sh` to write the Markdown table of contents. |
| `user_guide.css` | Extra styles for the guide page (loaded after `main.css` and `custom.css`). |

## Updating the guide

1. Edit `nerdss_user_guide.md`. It is GitHub-flavoured Markdown:
   - Use `##` for sections and `###` for subsections. These two levels make up the table of contents.
     Parameters are `####` (MOL file) or `#####` (INP file) headings, so each has its own link, for example
     `nerdss_user_guide.html#trajwrite`.
   - Put maths in `$...$` or `$$...$$` (TeX). The page renders it with MathJax.
   - Put images in `media/` and link them relatively: `![description](media/name.png)`.
   - Do not edit between `<!-- toc -->` and `<!-- tocstop -->`. `build.sh` rewrites that block.
2. Run the build (needs `pandoc` 3.x and `python3`):

   ```sh
   ./user_guide/build.sh
   ```

3. Commit the `.md` and the regenerated `.html` together.

## Where it came from

The Markdown was converted with pandoc from `NERDSS_USER_GUIDE_MASTER.docx`, last updated on 28 October 2026 on the
`nerdss-optimized-1.4` branch of the NERDSS repository. It was then cleaned up by hand. The text was kept as written,
but the formatting was changed:

- Word artefacts (bookmarks, underline spans, empty headings) were removed.
- Commands, file names and reaction rules are now code. En dashes that Word had put into commands
  (`–f`, `–r`, `–a`, `–j`) were turned back into hyphens.
- Sample blocks are now code blocks.
- The time-step equations are now displayed TeX.
- Sections that are not about input and output files are now their own top-level sections: "NERDSS parameters",
  "Use a coordinate file", "Rule-based reaction-rate considerations" and "Time-step selection".

The Word document is still the archival copy, and so is `files/NERDSS_USER_GUIDE.pdf`, an older 2024 edition.
From now on, make content changes here in Markdown.
