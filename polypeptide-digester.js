function cleave(sequence, enzymeCleavagePatterns, missed_cleavages, min_length, max_length) {
    // enzymeCleavagePatterns = [[cleavage site regex, (look-behind pattern or false)]]
    // javascript regular expression engine doesn't support 'look-behind'
    // patterns, which are necessary to match overlapping cleavage sites.
    // This script is using XRegExp module and xregexp-lookbehind2 plugin,
    // which overcomes this issue.
    // Each enzyme in enzymeRegExpArray dictionary contains at least one array with regex pattern.
    // Each regex array contains regex pattern stripped from look-behind expression
    // i.e. enzymeCleavagePatterns[0][0] and extracted look-behind pattern as a string
    // i.e. enzymeCleavagePatterns[0][1] or false when look-behind is not needed

    function sortNumber(a,b) {
        return a - b;
    }

    var cleavage_sites= [];
    for (pattern of enzymeCleavagePatterns){
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

    var missed_cleavages = Number(missed_cleavages);
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
    return nonRedundantPeptides;
}
