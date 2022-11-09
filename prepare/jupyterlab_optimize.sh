#%matplotlib widget
pip install jupyterlab-language-pack-zh-CN
pip install -U jedi-language-server
# conda install -c conda-forge texlab chktex
jupyter serverextension enable --py jupyterlab --sys-prefix
# pip install --upyterlab-git
# pip install --upgrade jupyterlab jupyterlab-git
pip install jupyterlab-github
jupyter labextension install @jupyterlab/toc
pip install jupyterlab_latex
jupyter labextension install @mflevine/jupyterlab_html
pip install jupyterlab-fasta
pip install jupyterlab-geojson
pip install jupyterlab-katex
pip install jupyterlab-mathjax3
pip install jupyterlab-vega2
pip install jupyterlab-vega3
pip install jupyter_bokeh
pip install jupyterlab-drawio
pip install jupyterlab_sql
jupyter serverextension enable jupyterlab_sql --py --sys-prefix
jupyter lab build
pip install lckr-jupyterlab-variableinspector
pip install "jupyterlab>=1.0" jupyterlab-dash==0.1.0a3
jupyter labextension install nbgather
pip install 'jupyterlab>=3.0.0,<4.0.0a0' jupyterlab-lsp
jupyter labextension install jupyterlab-spreadsheet
# https://zhuanlan.zhihu.com/p/437592449⬇
# pip install "jupyterlab-kite>=2.0.2"
pip install blackcellmagic
#%load_ext blackcellmagic
#%%black
pip install watermark
conda install -c conda-forge nodejs
#load_ext watermark
#%watermark
conda install jupyterlab -c conda-forge
