# ~/.bashrc: executed by bash(1) for non-login shells.
# see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
# for examples

# If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac

# don't put duplicate lines or lines starting with space in the history.
# See bash(1) for more options
HISTCONTROL=ignoreboth

# append to the history file, don't overwrite it
shopt -s histappend

# for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
HISTSIZE=1000
HISTFILESIZE=2000

# check the window size after each command and, if necessary,
# update the values of LINES and COLUMNS.
shopt -s checkwinsize

# If set, the pattern "**" used in a pathname expansion context will
# match all files and zero or more directories and subdirectories.
#shopt -s globstar

# make less more friendly for non-text input files, see lesspipe(1)
[ -x /usr/bin/lesspipe ] && eval "$(SHELL=/bin/sh lesspipe)"

# set variable identifying the chroot you work in (used in the prompt below)
if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

# set a fancy prompt (non-color, unless we know we "want" color)
case "$TERM" in
    xterm-color|*-256color) color_prompt=yes;;
esac

# uncomment for a colored prompt, if the terminal has the capability; turned
# off by default to not distract the user: the focus in a terminal window
# should be on the output of commands, not on the prompt
#force_color_prompt=yes

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
	# We have color support; assume it's compliant with Ecma-48
	# (ISO/IEC-6429). (Lack of such support is extremely rare, and such
	# a case would tend to support setf rather than setaf.)
	color_prompt=yes
    else
	color_prompt=
    fi
fi

if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
else
    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
fi
unset color_prompt force_color_prompt

# If this is an xterm set the title to user@host:dir
case "$TERM" in
xterm*|rxvt*)
    PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
    ;;
*)
    ;;
esac

# enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

# colored GCC warnings and errors
#export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'

# some more ls aliases
alias ll='ls -alF'
alias la='ls -A'
alias l='ls -CF'

# Add an "alert" alias for long running commands.  Use like so:
#   sleep 10; alert
alias alert='notify-send --urgency=low -i "$([ $? = 0 ] && echo terminal || echo error)" "$(history|tail -n1|sed -e '\''s/^\s*[0-9]\+\s*//;s/[;&|]\s*alert$//'\'')"'

# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.

if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi

# enable programmable completion features (you don't need to enable
# this, if it's already enabled in /etc/bash.bashrc and /etc/profile
# sources /etc/bash.bashrc).
if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fi

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/zhangxin/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/zhangxin/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/zhangxin/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/zhangxin/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


# >>> clash initialize >>>
#export https_proxy=http://127.0.0.1:7890 
#export http_proxy=http://127.0.0.1:7890
#export all_proxy=socks5://127.0.0.1:7890
# git config --global http.proxy 'http://127.0.0.1:7890'
# git config --global https.proxy 'https://127.0.0.1:7890'
# <<< clash initialize <<<

# >>> zhangxin alias >>>
alias zset.sh='bash ~/content/zset/zset.sh'
# <<< zhangxin alias <<<
# >>> zhangxin alias:2023-03-02-01-51-43 >>>
alias test.sh='bash /data/data/com.termux/files/home/content/zbash/test.sh'
alias up_.sh='bash /data/data/com.termux/files/home/content/zbash/up_.sh'
alias zbash.sh='bash /data/data/com.termux/files/home/content/zbash/zbash.sh'
alias zx0.sh='bash /data/data/com.termux/files/home/content/zbash/zx0.sh'
alias zx0X.sh='bash /data/data/com.termux/files/home/content/zbash/zx0X.sh'
alias zx1.sh='bash /data/data/com.termux/files/home/content/zbash/zx1.sh'
alias zx1X.sh='bash /data/data/com.termux/files/home/content/zbash/zx1X.sh'
alias zx2.sh='bash /data/data/com.termux/files/home/content/zbash/zx2.sh'
alias zx2X.sh='bash /data/data/com.termux/files/home/content/zbash/zx2X.sh'
alias zx8.sh='bash /data/data/com.termux/files/home/content/zbash/zx8.sh'
alias zx9.sh='bash /data/data/com.termux/files/home/content/zbash/zx9.sh'
# <<< zhangxin alias:2023-03-02-01-51-43 <<<
# >>> zhangxin alias:2023-03-02-01-55-06 >>>
alias gpush.sh='bash /data/data/com.termux/files/home/content/zbash/gpush.sh'
alias test.sh='bash /data/data/com.termux/files/home/content/zbash/test.sh'
alias up_.sh='bash /data/data/com.termux/files/home/content/zbash/up_.sh'
alias zbash.sh='bash /data/data/com.termux/files/home/content/zbash/zbash.sh'
alias zx0.sh='bash /data/data/com.termux/files/home/content/zbash/zx0.sh'
alias zx0X.sh='bash /data/data/com.termux/files/home/content/zbash/zx0X.sh'
alias zx1.sh='bash /data/data/com.termux/files/home/content/zbash/zx1.sh'
alias zx1X.sh='bash /data/data/com.termux/files/home/content/zbash/zx1X.sh'
alias zx2.sh='bash /data/data/com.termux/files/home/content/zbash/zx2.sh'
alias zx2X.sh='bash /data/data/com.termux/files/home/content/zbash/zx2X.sh'
alias zx8.sh='bash /data/data/com.termux/files/home/content/zbash/zx8.sh'
alias zx9.sh='bash /data/data/com.termux/files/home/content/zbash/zx9.sh'
# <<< zhangxin alias:2023-03-02-01-55-06 <<<
# >>> zhangxin alias:2023-03-02-01-56-45 >>>
alias gpush.sh='bash /data/data/com.termux/files/home/content/zbash/gpush.sh'
alias slog.sh='bash /data/data/com.termux/files/home/content/zbash/slog.sh'
alias test.sh='bash /data/data/com.termux/files/home/content/zbash/test.sh'
alias up_.sh='bash /data/data/com.termux/files/home/content/zbash/up_.sh'
alias zbash.sh='bash /data/data/com.termux/files/home/content/zbash/zbash.sh'
alias zx0.sh='bash /data/data/com.termux/files/home/content/zbash/zx0.sh'
alias zx0X.sh='bash /data/data/com.termux/files/home/content/zbash/zx0X.sh'
alias zx1.sh='bash /data/data/com.termux/files/home/content/zbash/zx1.sh'
alias zx1X.sh='bash /data/data/com.termux/files/home/content/zbash/zx1X.sh'
alias zx2.sh='bash /data/data/com.termux/files/home/content/zbash/zx2.sh'
alias zx2X.sh='bash /data/data/com.termux/files/home/content/zbash/zx2X.sh'
alias zx8.sh='bash /data/data/com.termux/files/home/content/zbash/zx8.sh'
alias zx9.sh='bash /data/data/com.termux/files/home/content/zbash/zx9.sh'
# <<< zhangxin alias:2023-03-02-01-56-45 <<<
# >>> zhangxin alias:2023-03-02-14-00-53 >>>
alias gpull.sh='bash /data/data/com.termux/files/home/content/zbash/gpull.sh'
alias gpush.sh='bash /data/data/com.termux/files/home/content/zbash/gpush.sh'
alias slog.sh='bash /data/data/com.termux/files/home/content/zbash/slog.sh'
alias test.sh='bash /data/data/com.termux/files/home/content/zbash/test.sh'
alias up_.sh='bash /data/data/com.termux/files/home/content/zbash/up_.sh'
alias zbash.sh='bash /data/data/com.termux/files/home/content/zbash/zbash.sh'
alias zx0.sh='bash /data/data/com.termux/files/home/content/zbash/zx0.sh'
alias zx0X.sh='bash /data/data/com.termux/files/home/content/zbash/zx0X.sh'
alias zx1.sh='bash /data/data/com.termux/files/home/content/zbash/zx1.sh'
alias zx1X.sh='bash /data/data/com.termux/files/home/content/zbash/zx1X.sh'
alias zx2.sh='bash /data/data/com.termux/files/home/content/zbash/zx2.sh'
alias zx2X.sh='bash /data/data/com.termux/files/home/content/zbash/zx2X.sh'
alias zx8.sh='bash /data/data/com.termux/files/home/content/zbash/zx8.sh'
alias zx9.sh='bash /data/data/com.termux/files/home/content/zbash/zx9.sh'
# <<< zhangxin alias:2023-03-02-14-00-53 <<<
# >>> zhangxin alias:2023-03-02-14-14-23 >>>
alias gpull.sh='bash /data/data/com.termux/files/home/content/zbash/gpull.sh'
alias gpush.sh='bash /data/data/com.termux/files/home/content/zbash/gpush.sh'
alias slog.sh='bash /data/data/com.termux/files/home/content/zbash/slog.sh'
alias test.sh='bash /data/data/com.termux/files/home/content/zbash/test.sh'
alias up_.sh='bash /data/data/com.termux/files/home/content/zbash/up_.sh'
alias zbash.sh='bash /data/data/com.termux/files/home/content/zbash/zbash.sh'
alias zx0.sh='bash /data/data/com.termux/files/home/content/zbash/zx0.sh'
alias zx0X.sh='bash /data/data/com.termux/files/home/content/zbash/zx0X.sh'
alias zx1.sh='bash /data/data/com.termux/files/home/content/zbash/zx1.sh'
alias zx1X.sh='bash /data/data/com.termux/files/home/content/zbash/zx1X.sh'
alias zx2.sh='bash /data/data/com.termux/files/home/content/zbash/zx2.sh'
alias zx2X.sh='bash /data/data/com.termux/files/home/content/zbash/zx2X.sh'
alias zx8.sh='bash /data/data/com.termux/files/home/content/zbash/zx8.sh'
alias zx9.sh='bash /data/data/com.termux/files/home/content/zbash/zx9.sh'
# <<< zhangxin alias:2023-03-02-14-14-23 <<<
# >>> zhangxin alias:2023-03-03-14-03-08 >>>
alias gpull.sh='bash /data/data/com.termux/files/home/content/zbash/gpull.sh'
alias gpush.sh='bash /data/data/com.termux/files/home/content/zbash/gpush.sh'
alias slog.sh='bash /data/data/com.termux/files/home/content/zbash/slog.sh'
alias test.sh='bash /data/data/com.termux/files/home/content/zbash/test.sh'
alias up_.sh='bash /data/data/com.termux/files/home/content/zbash/up_.sh'
alias zbash.sh='bash /data/data/com.termux/files/home/content/zbash/zbash.sh'
alias zx0.sh='bash /data/data/com.termux/files/home/content/zbash/zx0.sh'
alias zx0X.sh='bash /data/data/com.termux/files/home/content/zbash/zx0X.sh'
alias zx1.sh='bash /data/data/com.termux/files/home/content/zbash/zx1.sh'
alias zx1X.sh='bash /data/data/com.termux/files/home/content/zbash/zx1X.sh'
alias zx2.sh='bash /data/data/com.termux/files/home/content/zbash/zx2.sh'
alias zx2X.sh='bash /data/data/com.termux/files/home/content/zbash/zx2X.sh'
alias zx8.sh='bash /data/data/com.termux/files/home/content/zbash/zx8.sh'
alias zx9.sh='bash /data/data/com.termux/files/home/content/zbash/zx9.sh'
# <<< zhangxin alias:2023-03-03-14-03-08 <<<
