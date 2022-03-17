library('data.table')
library('lintr')
library('rex')

double_quotes_linter <- function(source_file) {
  content <- source_file$full_parsed_content
  str_idx <- which(content$token == "STR_CONST")
  squote_matches <- which(re_matches(
    content[str_idx, "text"],
    rex(start, double_quote, any_non_single_quotes, double_quote, end)
  ))

  lapply(
    squote_matches,
    function(id) {
      with(content[str_idx[id], ], {
        line <- source_file$file_lines[line1]
        col2 <- if (line1 == line2) col2 else nchar(line)
        Lint(
          filename = source_file$filename,
          line_number = line1,
          column_number = col1,
          type = "style",
          message = "Only use single-quotes.",
          line = line,
          ranges = list(c(col1, col2))
        )
      })
    }
  )
}
newDefaults = with_defaults(double_quotes_linter = double_quotes_linter,
                            line_length_linter(120),
                            assignment_linter = NULL,
                            closed_curly_linter = NULL,
                            object_name_linter = object_name_linter('camelCase'),
                            single_quotes_linter = NULL,
                            commented_code_linter = NULL)
lintsFound = lint_package(linters = newDefaults)
lfDt = unique(as.data.table(lintsFound), by = c('filename', 'line_number', 'message'))
lfDt[, lint_link := paste0('https://github.com/hugheylab/', repository, '/blob/', branch, '/', filename, '#L', line_number)]
lfDt[, line := trimws(line)]
setorder(lfDt, filename, line_number)

# %0D = \r and %0A = \n
# Needs different formatting for bash output
newlineEsc = ' %0D%0A'

lfDt[, format_line :=
       paste0('- ', filename, ' line ', line_number, ': ', message, ' (',
              lint_link, ')', newlineEsc, '```', newlineEsc, line, newlineEsc, '```')]
issueStr = paste0(lfDt$format_line, collapse = newlineEsc)

issueStr = gsub("'", "'\"'\"'", issueStr, fixed = TRUE)

s = sprintf("echo '::set-output name=style_text::%s'", issueStr)
