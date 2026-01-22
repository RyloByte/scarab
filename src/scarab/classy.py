__author__ = 'Ryan J McLaughlin'


class ScarabBase:
    """
    A base class for all SCARAB sub-commands. It requires shared properties
    """

    def __init__(self, subcmd_name) -> None:
        self.subcmd = subcmd_name
        self.executables = {}
        self.aln_file = ""
        self.seq_file = ""
        self.output_sep = ','
        return

    def get_info(self) -> str:
        info_string = ""
        return info_string

    def furnish_with_arguments(self, args) -> None:
        return
