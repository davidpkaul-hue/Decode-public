-- Convert <span style="color:green">…</span> into docx text runs with green color.
-- Pandoc's default docx writer ignores the CSS style attribute on Span, so we
-- emit raw OOXML for any Span whose style attribute contains "color: <name>".

local color_map = {
  green   = "008000",
  red     = "C00000",
  blue    = "0000FF",
  orange  = "ED7D31",
  yellow  = "BFBF00",
  purple  = "7030A0",
  gray    = "595959",
  black   = "000000",
}

local function hex_for(name)
  if not name then return nil end
  name = name:lower():gsub("%s+", "")
  if color_map[name] then return color_map[name] end
  -- #RRGGBB or RRGGBB
  local hex = name:match("^#?(%x%x%x%x%x%x)$")
  if hex then return hex:upper() end
  return nil
end

function Span(el)
  local style = el.attributes and el.attributes.style
  if not style then return nil end
  local color_value = style:match("color%s*:%s*([^;]+)")
  local hex = hex_for(color_value)
  if not hex then return nil end

  -- Render the span's inline content as OOXML and wrap each run with a color.
  local inner = pandoc.write(pandoc.Pandoc({pandoc.Plain(el.content)}), "docx")
  -- Fallback: if write returns binary docx, use stringify
  local text = pandoc.utils.stringify(el)

  -- Build a single colored run. Preserves whitespace via xml:space.
  local ooxml = string.format(
    '<w:r><w:rPr><w:color w:val="%s"/></w:rPr><w:t xml:space="preserve">%s</w:t></w:r>',
    hex,
    text:gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;")
  )
  return pandoc.RawInline("openxml", ooxml)
end
