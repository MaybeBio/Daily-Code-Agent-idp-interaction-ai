# Tab completion for noodle: press Tab after "noodle " for the commands, and after a
# command for the run names. setup.sh prints the line that loads this from your ~/.bashrc:
#     source <noodle-md folder>/noodle-completion.bash

_NOODLE_RUNS="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/runs"
_NOODLE_COMMANDS="new launch status analyze cancel export clean check help"

_noodle() {
    local word=${COMP_WORDS[COMP_CWORD]} command=${COMP_WORDS[1]} runs="" folder
    COMPREPLY=()
    if [ "$COMP_CWORD" -eq 1 ]; then
        COMPREPLY=($(compgen -W "$_NOODLE_COMMANDS" -- "$word"))
        return
    fi
    case "$command" in
        help)                                [ "$COMP_CWORD" -eq 2 ] || return
                                             COMPREPLY=($(compgen -W "$_NOODLE_COMMANDS" -- "$word"))
                                             return ;;
        launch|status|analyze|cancel|clean)  [ "$COMP_CWORD" -eq 2 ] || return ;;
        export)                              ;;  # any number of runs
        *)                                   return ;;
    esac
    for folder in "$_NOODLE_RUNS"/*/; do
        [ -e "$folder/run.yaml" ] && runs+=" $(basename "$folder")"
    done
    COMPREPLY=($(compgen -W "$runs" -- "$word"))
}

# also completes ./noodle and /any/path/noodle, and the alias noodle
complete -F _noodle noodle
