from __future__ import absolute_import
from __future__ import print_function
import re, fnmatch
from sys import stderr
from six.moves import zip

globalNuisances = re.compile(
    "(lumi|pdf_(qqbar|gg|qg)|QCDscale_(ggH|qqH|VH|ggH1in|ggH2in|VV)|UEPS|FakeRate|CMS_(eff|fake|trigger|scale|res)_([gemtjb]|met))"
)


from combine.Datacard import Datacard
from combine.NuisanceModifier import doEditNuisance


def isVetoed(name, vetoList):
    for pattern in vetoList:
        if not pattern:
            continue
        if re.match(pattern, name):
            return True
    return False


def isIncluded(name, includeList):
    if not len(includeList):
        return True
    for pattern in includeList:
        if not pattern:
            continue
        if re.match(pattern, name):
            return True
    return False


def addRateParam(lsyst, f, ret):

    if len(f) > 6 or len(f) < 5:
        raise RuntimeError(
            "Error, directives of type 'rateParam' should be of form .. name rateParam channel process initial value OR name rateParam channel process formula args"
        )
    if len(f) == 5:
        if ".root" in f[4] and ":" in f[4]:
            ty = 2
        else:
            ty = 0
        tmp_exp = [[lsyst, f[4], ty], ""]  # Case for free parameter with no range
    elif len(f) == 6:
        if "[" in f[-1] and "]" in f[-1]:
            tmp_exp = [[lsyst, f[4], 0], f[-1]]
        else:
            tmp_exp = [[lsyst, f[4], f[5], 1], ""]
    # check for malformed bin/process
    if f[2] not in ret.bins or f[3] not in ret.processes:
        raise RuntimeError(" No such channel/process '%s/%s', malformed line:\n   %s" % (f[2], f[3], " ".join(f)))
    if ("%sAND%s" % (f[2], f[3])) in list(ret.rateParams.keys()):
        ret.rateParams["%sAND%s" % (f[2], f[3])].append(tmp_exp)
    else:
        ret.rateParams["%sAND%s" % (f[2], f[3])] = [tmp_exp]
    ret.rateParamsOrder.add(lsyst)


def parseCard(file, bin, noJMax, allowNoSignal, allowNoBackground, evaluateEdits, nuisancesToExclude, stat, verbose):
    if type(file) == type("str"):
        raise RuntimeError(
            "You should pass as argument to parseCards a file object, stream or a list of lines, not a string"
        )
    ret = Datacard()

    # resetting these here to defaults, parseCard will fill them up
    ret.discretes = []
    ret.groups = {}

    #
    nbins = -1
    nprocesses = -1
    nuisances = -1
    binline = []
    processline = []
    sigline = []
    shapesUseBin = False
    lineNumber = None

    try:
        for lineNumber, l in enumerate(file):
            f = l.split()
            if len(f) < 1:
                continue
            if f[0] == "imax":
                nbins = int(f[1]) if f[1] != "*" else -1
            if f[0] == "jmax":
                nprocesses = int(f[1]) + 1 if f[1] != "*" else -1
            if f[0] == "kmax":
                nuisances = int(f[1]) if f[1] != "*" else -1
            if f[0] == "shapes":
                if not bin:
                    raise RuntimeError("Can use shapes only with binary output mode")
                if len(f) < 4:
                    raise RuntimeError("Malformed shapes line")
                if f[2] not in ret.shapeMap:
                    ret.shapeMap[f[2]] = {}
                if f[1] in ret.shapeMap[f[2]]:
                    raise RuntimeError("Duplicate definition for process '%s', channel '%s'" % (f[1], f[2]))
                ret.shapeMap[f[2]][f[1]] = f[3:]
                if "$CHANNEL" in l:
                    shapesUseBin = True
                if f[2] != "*":
                    shapesUseBin = True
            if f[0] == "Observation" or f[0] == "observation":
                ret.obs = [float(x) for x in f[1:]]
                if nbins == -1:
                    nbins = len(ret.obs)
                if len(ret.obs) != nbins:
                    raise RuntimeError("Found %d observations but %d bins have been declared" % (len(ret.obs), nbins))
                if binline != []:
                    if len(binline) != len(ret.obs):
                        raise RuntimeError(
                            "Found %d bins (%s) but %d bins have been declared" % (len(ret.bins), ret.bins, nbins)
                        )
                    ret.bins = binline
                    ret.obs = dict([(b, ret.obs[i]) for i, b in enumerate(ret.bins)])
                    binline = []
            if f[0] == "bin":
                binline = []
                for b in f[1:]:
                    if re.match("[0-9]+", b):
                        raise RuntimeError("Error: Bin %(b)s starts with a digit!" % locals())
                    binline.append(b)
            if f[0] == "process":
                if processline == []:  # first line contains names
                    processline = f[1:]
                    if len(binline) != len(processline):
                        raise RuntimeError("'bin' line has a different length than 'process' line.")
                    continue
                sigline = f[1:]  # second line contains ids
                if re.match("-?[0-9]+", processline[0]) and not re.match("-?[0-9]+", sigline[0]):
                    (processline, sigline) = (sigline, processline)
                if len(sigline) != len(processline):
                    raise RuntimeError("'bin' line has a different length than 'process' line.")
                hadBins = len(ret.bins) > 0
                for i, b in enumerate(binline):
                    p = processline[i]
                    s = int(sigline[i]) <= 0  # <=0 for signals, >0 for backgrounds
                    ret.keyline.append((b, processline[i], s))
                    if hadBins:
                        if b not in ret.bins:
                            raise RuntimeError("Bin %s not among the declared bins %s" % (b, ret.bins))
                    else:
                        if b not in ret.bins:
                            ret.bins.append(b)
                    if p not in ret.processes:
                        ret.processes.append(p)
                if nprocesses == -1:
                    nprocesses = len(ret.processes)
                if nbins == -1:
                    nbins = len(ret.bins)
                if not noJMax:
                    if nprocesses != len(ret.processes):
                        raise RuntimeError(
                            "Found %d processes (%s), declared jmax = %d"
                            % (len(ret.processes), ret.processes, nprocesses)
                        )
                if nbins != len(ret.bins):
                    raise RuntimeError("Found %d bins (%s), declared imax = %d" % (len(ret.bins), ret.bins, nbins))
                ret.exp = dict([(b, {}) for b in ret.bins])
                ret.isSignal = dict([(p, None) for p in ret.processes])
                if ret.obs != [] and type(ret.obs) == list:  # still as list, must change into map with bin names
                    ret.obs = dict([(b, ret.obs[i]) for i, b in enumerate(ret.bins)])
                for (b, p, s) in ret.keyline:
                    if ret.isSignal[p] == None:
                        ret.isSignal[p] = s
                    elif ret.isSignal[p] != s:
                        raise RuntimeError(
                            "Process %s is declared as signal in some bin and as background in some other bin" % p
                        )
                ret.signals = [p for p, s in ret.isSignal.items() if s == True]
                if len(ret.signals) == 0 and not allowNoSignal:
                    raise RuntimeError("You must have at least one signal process (id <= 0)")
            if f[0] == "rate":
                if processline == []:
                    raise RuntimeError("Missing line with process names before rate line")
                if sigline == []:
                    raise RuntimeError("Missing line with process id before rate line")
                if len(f[1:]) != len(ret.keyline):
                    raise RuntimeError(
                        "Malformed rate line: length %d, while bins and process lines have length %d"
                        % (len(f[1:]), len(ret.keyline))
                    )
                for (b, p, s), r in zip(ret.keyline, f[1:]):
                    ret.exp[b][p] = float(r)
                break  # rate is the last line before nuisances
        # parse nuisances
        for lineNumber, l in enumerate(file):
            if l.startswith("--"):
                continue
            l = re.sub("\\s*#.*", "", l)
            l = re.sub("(?<=\\s)-+(\\s|$)", " 0\\1", l)
            f = l.split()
            if len(f) <= 1:
                continue
            nofloat = False
            lsyst = f[0]
            pdf = f[1]
            args = []
            numbers = f[2:]
            if lsyst.endswith("[nofloat]"):
                lsyst = lsyst.replace("[nofloat]", "")
                nofloat = True
            if nuisancesToExclude and isVetoed(lsyst, nuisancesToExclude):
                if verbose > 0:
                    stderr.write(
                        "Excluding nuisance %s selected by a veto pattern among %s\n" % (lsyst, nuisancesToExclude)
                    )
                if nuisances != -1:
                    nuisances -= 1
                continue
            if re.match("[0-9]+", lsyst):
                lsyst = "theta" + lsyst
            if pdf == "lnN" or pdf == "lnU" or pdf == "gmM" or pdf == "trG" or pdf.startswith("shape"):
                pass  # nothing special to do
            elif pdf == "gmN":
                args = [int(f[2])]
                numbers = f[3:]
            elif pdf == "unif":
                args = [float(f[2]), float(f[3])]
                numbers = f[4:]
            elif pdf == "dFD" or pdf == "dFD2":
                args = [float(f[2])]
                numbers = f[3:]
            elif pdf == "param":
                # for parametric uncertainties, there's no line to account per bin/process effects
                # just assume everything else is an argument and move on
                args = f[2:]
                if len(args) <= 1:
                    raise RuntimeError(
                        "Uncertainties of type 'param' must have at least two arguments (mean and sigma)"
                    )
                ret.systs.append([lsyst, nofloat, pdf, args, []])
                continue
            elif pdf == "flatParam":
                ret.flatParamNuisances[lsyst] = True
                # for flat parametric uncertainties, code already does the right thing as long as they are non-constant RooRealVars linked to the model
                continue
            elif pdf == "extArg":
                # look for additional parameters in workspaces
                ret.extArgs[lsyst] = f[:]
                continue
            elif pdf == "rateParam":
                if ("*" in f[3]) or ("*" in f[2]):  # all channels/processes
                    for c in ret.processes:
                        for b in ret.bins:
                            if not fnmatch.fnmatch(c, f[3]):
                                continue
                            if not fnmatch.fnmatch(b, f[2]):
                                continue
                            f_tmp = f[:]
                            f_tmp[2] = b
                            f_tmp[3] = c
                            addRateParam(lsyst, f_tmp, ret)
                else:
                    addRateParam(lsyst, f, ret)
                continue
            elif pdf == "discrete":
                args = f[2:]
                ret.discretes.append(lsyst)
                continue
            elif pdf == "edit":
                if nuisances != -1:
                    nuisances = -1
                if evaluateEdits:
                    if verbose > 1:
                        print("Before edit: \n\t%s\n" % ("\n\t".join([str(x) for x in ret.systs])))
                    if verbose > 1:
                        print("Edit command: %s\n" % numbers)
                    doEditNuisance(ret, numbers[0], numbers[1:])
                    if verbose > 1:
                        print("After edit: \n\t%s\n" % ("\n\t".join([str(x) for x in ret.systs])))
                else:
                    if numbers[0] in ["changepdf", "freeze"]:
                        ret.nuisanceEditLines.append([numbers[0], numbers[1:]])
                    else:
                        ret.nuisanceEditLines.append([numbers[0], numbers[1], numbers[2], numbers[3:]])
                continue
            elif pdf == "group":
                # This is not really a pdf type, but a way to be able to name groups of nuisances together
                groupName = lsyst
                groupNuisances = numbers

                if not groupNuisances:
                    raise RuntimeError("Syntax error for group '%s': empty line after 'group'." % groupName)

                defToks = ("=", "+=")
                defTok = groupNuisances.pop(0)
                if defTok not in defToks:
                    raise RuntimeError(
                        "Syntax error for group '%s': first thing after 'group' is not '[+]=' but '%s'."
                        % (groupName, defTok)
                    )

                if groupName not in ret.groups:
                    if defTok == "=":
                        ret.groups[groupName] = set(groupNuisances)
                    else:
                        raise RuntimeError("Cannot append to group '%s' as it was not yet defined." % groupName)
                else:
                    if defTok == "+=":
                        ret.groups[groupName].update(set(groupNuisances))
                    else:
                        raise RuntimeError(
                            "Will not redefine group '%s'. It previously contained '%s' and you now wanted it to contain '%s'."
                            % (groupName, ret.groups[groupName], groupNuisances)
                        )

                continue
            elif pdf == "autoMCStats":
                if len(f) > 5:
                    raise RuntimeError(
                        "Syntax for autoMCStats should be 'channel autoMCStats threshold [include-signal = 0] [hist-mode = 0]"
                    )
                statThreshold = float(f[2])
                statIncludeSig = bool(int(f[3])) if len(f) >= 4 else False
                statHistMode = int(f[4]) if len(f) >= 5 else 1
                statFlags = (statThreshold, statIncludeSig, statHistMode)
                if "*" in lsyst:
                    for b in ret.bins:
                        if not fnmatch.fnmatch(b, lsyst):
                            continue
                        ret.binParFlags[b] = statFlags
                else:
                    if lsyst not in ret.bins:
                        raise RuntimeError(" No such channel '%s', malformed line:\n   %s" % (lsyst, " ".join(f)))
                    ret.binParFlags[lsyst] = statFlags
                continue
            else:
                raise RuntimeError("Unsupported pdf %s" % pdf)
            if len(numbers) < len(ret.keyline):
                raise RuntimeError(
                    "Malformed systematics line %s of length %d: while bins and process lines have length %d"
                    % (lsyst, len(numbers), len(ret.keyline))
                )
            errline = dict([(b, {}) for b in ret.bins])
            nonNullEntries = 0
            for (b, p, s), r in zip(ret.keyline, numbers):
                if "/" in r:  # "number/number"
                    if (pdf not in ["lnN", "lnU"]) and ("?" not in pdf):
                        raise RuntimeError("Asymmetric errors are allowed only for Log-normals")
                    errline[b][p] = [float(x) for x in r.split("/")]
                    for v in errline[b][p]:
                        if v <= 0.00:
                            raise ValueError(
                                'Found "%s" in the nuisances affecting %s for %s. This would lead to NANs later on, so please fix it.'
                                % (r, p, b)
                            )
                else:
                    errline[b][p] = float(r)
                    # values of 0.0 are treated as 1.0; scrap negative values.
                    if pdf not in ["trG", "dFD", "dFD2"] and errline[b][p] < 0:
                        raise ValueError(
                            'Found "%s" in the nuisances affecting %s in %s. This would lead to NANs later on, so please fix it.'
                            % (r, p, b)
                        )
                # set the rate to epsilon for backgrounds with zero observed sideband events.
                if pdf == "gmN" and ret.exp[b][p] == 0 and float(r) != 0:
                    ret.exp[b][p] = 1e-6
            ret.systs.append([lsyst, nofloat, pdf, args, errline])
    except Exception as ex:
        if lineNumber != None:
            msg = "Error reading line %d" % (lineNumber + 1)
            if hasattr(file, "name"):
                msg += " of file " + file.name

            msg += ": " + ex.args[0]
            ex.args = (msg,) + ex.args[1:]

        raise

    # check if there are bins with no rate
    for b in ret.bins:
        np_bin = sum([(ret.exp[b][p] != 0) for (b1, p, s) in ret.keyline if b1 == b])
        ns_bin = sum([(ret.exp[b][p] != 0) for (b1, p, s) in ret.keyline if b1 == b and s == True])
        nb_bin = sum([(ret.exp[b][p] != 0) for (b1, p, s) in ret.keyline if b1 == b and s != True])
        if np_bin == 0:
            raise RuntimeError("Bin %s has no processes contributing to it" % b)
        if ns_bin == 0 and not allowNoSignal:
            stderr.write("Warning: Bin %s has no signal processes contributing to it\n" % b)
        if nb_bin == 0 and not allowNoBackground:
            raise RuntimeError("Bin %s has no background processes contributing to it" % b)
    # cleanup systematics that have no effect to avoid zero derivatives
    syst2 = []
    for lsyst, nofloat, pdf, args, errline in ret.systs:
        nonNullEntries = 0
        if pdf == "param" or pdf == "discrete" or pdf == "rateParam":  # this doesn't have an errline
            syst2.append((lsyst, nofloat, pdf, args, errline))
            continue
        for (b, p, s) in ret.keyline:
            r = errline[b][p]
            nullEffect = r == 0.0 or (pdf == "lnN" and r == 1.0)
            if not nullEffect and ret.exp[b][p] != 0:
                nonNullEntries += 1  # is this a zero background?
        if nonNullEntries != 0:
            syst2.append((lsyst, nofloat, pdf, args, errline))
        elif nuisances != -1:
            nuisances -= 1  # remove from count of nuisances, since qe skipped it
    ret.systs = syst2
    # remove them if options.stat asks so
    if stat:
        nuisances = 0
        ret.systs = []
    # check number of nuisances
    if nuisances == -1:
        nuisances = len(ret.systs)
    elif len(ret.systs) != nuisances:
        raise RuntimeError("Found %d systematics, expected %d" % (len(ret.systs), nuisances))
    # set boolean to know about shape
    ret.hasShapes = len(ret.shapeMap) > 0
    # return result
    return ret


def FloatToString(inputValue):
    return ("%.10f" % inputValue).rstrip("0").rstrip(".")


def FloatToStringScientific(inputValue, etype="g"):
    s = FloatToString(inputValue)
    q = s.replace(".", "").lstrip("0")
    nq = max(len(q) - 1, 4)
    strcmd = "%%%c%i%c" % (".", nq, etype)
    return strcmd % inputValue
