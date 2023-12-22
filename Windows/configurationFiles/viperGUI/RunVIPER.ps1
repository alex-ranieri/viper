wsl --cd $HOME -d VIPER -e bash -c @'
 exec bash --rcfile <(echo 'source ~/.bashrc && python ~/viperGUI/viperGUI.py && exit')
'@