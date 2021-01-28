#git init
git add .
#git commit -m "first commit"
VERSION=$(date +'%H%M/%m%d%Y')
git commit -m "$VERSION"
#git branch -M main
#git remote add origin https://github.com/mourelatos-lab/TERA-Seq_manuscript.git
git push -u origin main
