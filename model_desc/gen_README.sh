# for now I am using TeXify on GitHub
python3 -m readme2tex  --svgdir model_desc/svgs/ --output README.md README.tex.md
mv model_desc/svgs .
rmdir model_desc


