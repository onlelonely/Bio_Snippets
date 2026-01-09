#!/bin/bash
# ---------------------------------------------
# Title: AI-Scientist-v2
# Description: From: Source/1. Atlas/Reference/Repo/AI-Scientist-v2.md
# ---------------------------------------------

Git clone https://github.com/SakanaAI/AI-Scientist-v2.git 
cd AI-Scientist-v2
pip install -R requirements.txt
export OPENAI_API_KEY='your_api_key_here'
python -m ai_scientist --task ****** --prompt "Rule 30"