#!/bin/bash
# ---------------------------------------------
# Title: CLI code
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Computational Tools/CLI code.md
# ---------------------------------------------

# 安裝nvm已安裝Node js > v18
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.3/install.sh | bash
# 記得restart
# 這邊裝v22
nvm install 22
nvm use 22 | npm install -g @anthropic-ai/claude-code
npm install -g @google/gemini-cli
npm install -g @openai/codex

#開啟Claude code
claude