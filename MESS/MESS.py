from util import *

import logging
LOGGER = logging.getLogger(__name__)

class MESS(object):
    """ A Massive eco-evolutionary synthesis simulation object.

    """

    def __init__(self, name, quiet=False, **kwargs):
        if not name:
            raise MESSError(REQUIRE_NAME)

        ## Do some checking here to make sure the name doesn't have
        ## special characters, spaces, or path delimiters. Allow _ and -.
        ## This will raise an error immediately if there are bad chars in name.
        self._check_name(name)


    ## Test assembly name is valid and raise if it contains any special characters
    def _check_name(self, name):
        invalid_chars = string.punctuation.replace("_", "")\
                                          .replace("-", "")+ " "
        if any(char in invalid_chars for char in name):
            raise MESSError(BAD_MESS_NAME.format(name))


    def write_params(self, outfile=None, force=False):
        """ Write out the parameters of this model to a file properly
        formatted as input for `mess -p <params.txt>`. A good and
        simple way to share/archive parameter settings for simulations.
        This is also the function that's used by __main__ to
        generate default params.txt files for `mess -n`
        """
        if outfile is None:
            outfile = "params-"+self.name+".txt"

        ## Test if params file already exists?
        ## If not forcing, test for file and bail out if it exists
        if not force:
            if os.path.isfile(outfile):
                raise IPyradWarningExit(PARAMS_EXISTS.format(outfile))

        with open(outfile, 'w') as paramsfile:
            ## Write the header. Format to 80 columns
            header = "------- mess params file (v.{})".format(MESS.__version__)
            header += ("-"*(80-len(header)))
            paramsfile.write(header)

            ## Whip through the current paramsdict and write out the current
            ## param value, the ordered dict index number. Also,
            ## get the short description from paramsinfo. Make it look pretty,
            ## pad nicely if at all possible.
            for key, val in self.paramsdict.iteritems():
                ## If multiple elements, write them out comma separated
                if isinstance(val, list) or isinstance(val, tuple):
                    paramvalue = ", ".join([str(i) for i in val])
                else:
                    paramvalue = str(val)

                ## skip deprecated params
                if key in ["edit_cutsites", "trim_overhang"]:
                    continue

                padding = (" "*(30-len(paramvalue)))
                paramkey = self.paramsdict.keys().index(key)
                paramindex = " ## [{}] ".format(paramkey)
                LOGGER.debug(key, val, paramindex)
                name = "[{}]: ".format(paramname(paramkey))
                description = paraminfo(paramkey, short=True)
                paramsfile.write("\n" + paramvalue + padding + \
                                        paramindex + name + description)

BAD_MESS_NAME = """\
    No spaces or special characters of any kind are allowed in the simulation 
    name. Special characters include all punctuation except dash '-' and 
    underscore '_'. A good practice is to replace spaces with underscores '_'.
    An example of a good simulation name is: hawaiian_arthropods 
    
    Here's what you put:
    {}
    """

REQUIRE_NAME = """\
    Simulation scenario name _must_ be set. This is the first parameter in the
    params.txt file, and will be used as a prefix for output files. It should be a 
    short string with no special characters, i.e., not a path (no \"/\" characters).
    If you need a suggestion, name it after the location you're working on.
    """

