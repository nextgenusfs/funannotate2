import sys
import os
import argparse
import shutil
from .__version__ import __version__
import textwrap as _textwrap


END_FORMATTING = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
SUPPRESS = "==SUPPRESS=="
OPTIONAL = "?"
ZERO_OR_MORE = "*"
ONE_OR_MORE = "+"
PARSER = "A..."
REMAINDER = "..."
_UNRECOGNIZED_ARGS_ATTR = "_unrecognized_args"


class MyParser(argparse.ArgumentParser):
    """
    This subclass of ArgumentParser changes the error messages, such that if a command is run with
    no other arguments, it will display the help text. If there is a different error, it will give
    the normal response (usage and error).
    """

    def error(self, message):
        if len(sys.argv) == 2:  # if a command was given but nothing else
            self.print_help(file=sys.stderr)
            sys.exit(2)
        else:
            super().error(message)


class MyHelpFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ["COLUMNS"] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if action.default != argparse.SUPPRESS and action.default is not None:
            if "default" not in help_text.lower():
                help_text += " (default: {})".format(action.default)
            elif "default: DEFAULT" in help_text:
                help_text = help_text.replace(
                    "default: DEFAULT", "default: {}".format(action.default)
                )
        return help_text

    def _split_lines(self, text, width):
        """
        Override this method to add special behaviour for help texts that start with:
          'R|' - loop text one option per line
        """
        return argparse.HelpFormatter._split_lines(self, text, width)

    def _format_action(self, action):
        result = super(MyHelpFormatter, self)._format_action(action)
        if isinstance(action, argparse._SubParsersAction):
            return "%*s%s" % (self._current_indent, "", result.lstrip())
        return result

    def _format_action_invocation(self, action):
        if isinstance(action, argparse._SubParsersAction):
            return ""
        return super(MyHelpFormatter, self)._format_action_invocation(action)

    def _iter_indented_subactions(self, action):
        if isinstance(action, argparse._SubParsersAction):
            try:
                get_subactions = action._get_subactions
            except AttributeError:
                pass
            else:
                for subaction in get_subactions():
                    yield subaction
        else:
            for subaction in super(MyHelpFormatter, self)._iter_indented_subactions(
                action
            ):
                yield subaction

    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(" ", text).strip()
        paragraphs = text.split("|n")
        multiline_text = ""
        for paragraph in paragraphs:
            formatted_paragraph = (
                _textwrap.fill(
                    paragraph, width, initial_indent=indent, subsequent_indent=indent
                )
                + "\n"
            )
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

    def _format_args(self, action, default_metavar):
        get_metavar = self._metavar_formatter(action, default_metavar)
        if action.nargs is None:
            result = "%s" % get_metavar(1)
        elif action.nargs == OPTIONAL:
            result = "[%s]" % get_metavar(1)
        elif action.nargs == ZERO_OR_MORE:
            metavar = get_metavar(1)
            if len(metavar) == 2:
                result = "[%s [%s ...]]" % metavar
            else:
                result = "[%s ...]" % metavar
        elif action.nargs == ONE_OR_MORE:
            result = "%s" % get_metavar(1)
        elif action.nargs == REMAINDER:
            result = "..."
        elif action.nargs == PARSER:
            result = "%s ..." % get_metavar(1)
        elif action.nargs == SUPPRESS:
            result = ""
        else:
            try:
                formats = ["%s" for _ in range(action.nargs)]
            except TypeError:
                raise ValueError("invalid nargs value") from None
            result = " ".join(formats) % get_metavar(action.nargs)
        return result
