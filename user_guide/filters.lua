-- Pandoc Lua filter for the NERDSS user guide (used by build.sh).
--  * the first level-1 heading becomes the document title (and is removed
--    from the body), so the sections start at level 2;
--  * tables are wrapped in <div class="table-wrapper"> so that wide tables
--    scroll sideways inside the page instead of widening it;
--  * GitHub-style alerts ("> [!WARNING]") become styled callouts.

local function is_html()
  return FORMAT:match('html') ~= nil
end

function Pandoc(doc)
  local first = doc.blocks[1]
  if first and first.t == 'Header' and first.level == 1 then
    if not doc.meta.title then
      doc.meta.title = pandoc.MetaInlines(first.content)
    end
    doc.blocks:remove(1)
  end
  return doc
end

function Table(tbl)
  if not is_html() then return nil end
  return pandoc.Div({ tbl }, { class = 'table-wrapper' })
end

local alert_titles = {
  NOTE = 'Note', TIP = 'Tip', IMPORTANT = 'Important',
  WARNING = 'Warning', CAUTION = 'Caution',
}

function BlockQuote(bq)
  if not is_html() then return nil end
  local first = bq.content[1]
  if not (first and first.t == 'Para' and first.content[1]
          and first.content[1].t == 'Str') then
    return nil
  end
  local kind = first.content[1].text:match('^%[!(%u+)%]$')
  if not kind or not alert_titles[kind] then return nil end
  first.content:remove(1)
  while first.content[1] and (first.content[1].t == 'SoftBreak'
        or first.content[1].t == 'LineBreak' or first.content[1].t == 'Space') do
    first.content:remove(1)
  end
  local blocks = pandoc.List(bq.content)
  if #first.content == 0 then blocks:remove(1) end
  return pandoc.Div(blocks, { class = 'guide-alert guide-alert-' .. kind:lower(),
                              role = 'note' })
end
