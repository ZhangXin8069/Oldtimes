#Begin.

#Install jupyterlab.
cd ~
conda install -c conda-forge jupyterlab
jupyter notebook --generate-config
echo 'import webbrowser
webbrowser.register('edge',None,webbrowser.GenericBrowser('/usr/bin/edge'))
c.NotebookApp.browser = 'edge'
c.NotebookApp.use_redirect_file = Fals
c.NotebookApp.allow_root = True' >> .jupyter/jupyter_notebook_config.py
echo "export BROWSER='/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe'" >> .bashrc
source .bashrc
jupyter lab
#jupyter-lab --allow-root --no-browser
#Instead, copy the link directly from the command line.

#End
