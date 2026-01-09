#!/bin/bash
# ---------------------------------------------
# Title: open_deep_research
# Description: From: Source/1. Atlas/Reference/Repo/open_deep_research.md
# ---------------------------------------------

Git clone https://github.com/langchain-ai/open_deep_research.git 
cd open_deep_research
pip install -R requirements.txt

# TAVILY_API_KEY="your_tavily_api_key"
# OPENAI_API_KEY="your_openai_api_key"
# PENAI_MODEL_NAME="gpt-4-turbo"

python -m deep_research --topic "The role of the gene [Your Gene of Interest] in [Your Disease of Interest]"