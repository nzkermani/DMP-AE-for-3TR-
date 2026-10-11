## Add to the project's .Rprofile (usethis::edit_r_profile(scope = "project"))
## so every RStudio session in this project is reachable by Claude Code's
## r-longomics MCP server. Claude Code connects to the session whose working
## directory matches the project it was launched in.
if (interactive() && requireNamespace("mcptools", quietly = TRUE)) {
  try(mcptools::mcp_session(), silent = TRUE)
}
