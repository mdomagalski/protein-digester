function cleave(sequence, enzyme, missed_cleavages, min_length, max_length) {
    //{enzyme_name: [cleavage site regex, (look-behind pattern or false)]}
    // javascript regular expression engine doesn't support 'look-behind'
    // patterns, which are necessary to match overlapping cleavage sites.
    // This script is using XRegExp module and xregexp-lookbehind2 plugin,
    // which overcomes this issue.
    // Each enzyme in expasy_rules dictionary contains at least one array with regex pattern.
    // Each regex array contains regex pattern stripped from look-behind expression
    // i.e. expasy_rules[enzyme][0][0] and extracted look-behind pattern as a string
    // i.e. expasy_rules[enzyme][0][1] or false when look-behind is not needed

    expasy_rules = {
        'arg-c':         [[/R/g, false]],
        'asp-n':         [[/\w(?=D)/g, false]],
        'bnps-skatole':  [[/W/g, false]],
        'caspase 1':     [[/D(?=[^PEDQKR])/g, "(?<=[FWYL]\w[HAT])"]],
        'caspase 2':     [[/D(?=[^PEDQKR])/g, "(?<=DVA)"]],
        'caspase 3':     [[/D(?=[^PEDQKR])/g, "(?<=DMQ)"]],
        'caspase 4':     [[/D(?=[^PEDQKR])/g, "(?<=LEV)"]],
        'caspase 5':     [[/D/g, "(?<=[LW]EH)"]],
        'caspase 6':     [[/D(?=[^PEDQKR])/g, "(?<=VE[HI])"]],
        'caspase 7':     [[/D(?=[^PEDQKR])/g, "(?<=DEV)"]],
        'caspase 8':     [[/D(?=[^PEDQKR])/g, "(?<=[IL]ET)"]],
        'caspase 9':     [[/D/g, "(?<=LEH)"]],
        'caspase 10':    [[/D/g, "(?<=IEA)"]],
        'chymotrypsin high specificity': [[/([FY](?=[^P]))|(W(?=[^MP]))/g, false]],
        'chymotrypsin low specificity': [[/([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))/g, false]],
        'clostripain':   [[/R/g, false]],
        'cnbr':          [[/M/g, false]],
        'enterokinase':  [[/K/g, "(?<=[DE]{3})"]],
        'factor xa':     [[/R/g, "(?<=[AFGILTVM][DE]G)"]],
        'formic acid':   [[/D/g, false]],
        'glutamyl endopeptidase': [[/E/g, false]],
        'granzyme b':    [[/D/g, "(?<=IEP)"]],
        'hydroxylamine': [[/N(?=G)/g, false]],
        'iodosobenzoic acid': [[/W/g, false]],
        'lysc':          [[/K/g, false]],
        'ntcb':          [[/\w(?=C)/g, false]],
        'pepsin ph2.0':  [[/([^R](?=[FLWY][^P]))|([FLWY](?=\w[^P]))/g, "(?<=[^HKR][^P])"]],
        'pepsin ph1.3':  [[/([^R](?=[FL][^P]))/g, "(?<=[^HKR][^P])"], [/([FL](?=\w[^P]))/g, "(?<=[^HKR][^P])"]],
        'proline endopeptidase': [[/P(?=[^P])/g, "(?<=[HKR])"]],
        'proteinase k':  [[/[AEFILTVWY]/g, false]],
        'staphylococcal peptidase i': [[/E/g, "(?<=[^E])"]],
        'thermolysin':   [[/[^DE](?=[AFILMV])/g, false]],
        'thrombin':      [[/(R(?=G))/, "(?<=G)"], [/(R(?=[^DE][^DE]))/g, "(?<=[AFGILTVM][AFGILTVWA]P)"]],
        'trypsin':       [[/([KR](?=[^P]))/g, false], [/K(?=P)/g, "(?<=W)"],[/R(?=P)/g, "(?<=M)"]]
    }

    function sortNumber(a,b) {
        return a - b;
    }

    if(!expasy_rules[enzyme]){
        throw "Not recognized enzyme: "+enzyme;
    }

    var cleavage_sites= [];
    for (pattern of expasy_rules[enzyme]){
        var regexpPattern = pattern[0];
        var lookbehindPattern = pattern[1];
        if (lookbehindPattern){
            var match=XRegExp.matchAllLb(sequence, lookbehindPattern, regexpPattern);
            for (var position of match){
                cleavage_sites.push(position);
            }
        }else{
            while(match=regexpPattern.exec(sequence)) {
                cleavage_sites.push(match.index + match[0].length);
            }
        }
    }
    cleavage_sites.sort(sortNumber);
    cleavage_sites.push(sequence.length);
    var cleavage_sites = cleavage_sites.filter(function(item, pos) {
        return cleavage_sites.indexOf(item) == Math.abs(pos);
    });

    var cleavage_sites_iterator = [0];
    var peptides = [];
    var cleavage = 1;
    var z =0;
    var range = Array.apply(null, Array(missed_cleavages+2)).map(function (_, i) {return i;});
    for(var i in cleavage_sites){
        cleavage_sites_iterator.push(cleavage_sites[i]);
        z+=1;
        if (cleavage < missed_cleavages+2){
            cleavage += 1;
        }
        for (var idx in range.slice(0,cleavage-1)){
            var seq = sequence.substring(cleavage_sites_iterator[cleavage_sites_iterator.length-cleavage+Number(idx)],cleavage_sites[i]);
            if ((min_length==null || seq.length >= min_length)&&(max_length==null ||seq.length <= max_length)){
                peptides.push(seq);
            }
        }
    }
    var nonRedundantPeptides = peptides.filter(function(item, pos) {
        return peptides.indexOf(item) == Math.abs(pos);
    });
    //console.log("sites: "+Number(z-1));
    //console.log("redundant pep count: "+peptides.length);
    //console.log("non-redundant pep count: "+nonRedundantPeptides.length);
    //console.log("peptides: "+peptides);
    return nonRedundantPeptides;
}
