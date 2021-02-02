# for now I am using TeXify on GitHub
python3 -m readme2tex  --output README.md README.tex.md
rm -rf ../svgs
mv svgs ..


