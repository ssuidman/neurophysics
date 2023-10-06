using CSV
using DataFrames

# include("store_data.jl")
# # Run this to store the variables in a data file called data.csv
# file_name = "data" #without .csv extension 
# variables = Dict("previn"=>previn,"pfminneg"=>pfminneg,"actualdiseases"=>actualdiseases,"pfmin"=>pfmin)
# store_variables(variables,file_name)

# # Run this to retrieve the variables back from the file data.csv
# file_name = "data" #without .csv extension 
# variables = get_variables(file_name) # is a dictionary with the name of the variables and the vectors of matrix that are stored in it
# previn,pfminneg,actualdiseases,pfmin = variables["previn"], variables["pfminneg"], variables["actualdiseases"], variables["pfmin"]' # to retrieve the actual variables

function store_variables(variables,file_name)
    df = DataFrame([])
    for var_name in keys(variables)
        if ndims(variables[var_name])==1
            df[!,var_name] = variables[var_name]
        elseif ndims(variables[var_name])==2
            for i in 1:size(variables[var_name], 1)
                df[!, "$(var_name) $i"] = variables[var_name][i, :]
            end
        end
    end
    CSV.write("/data/$file_name.csv",df) #without .csv extension 
end

function get_variables(file_name)
    df = DataFrame(CSV.File("$file_name.csv"))
    columns_df = names(df)
    variable_names = Dict{String, Vector{String}}()
    for (i,column) in enumerate(columns_df)
        var_name = split(column)[1]
        if !haskey(variable_names,var_name)
            variable_names[var_name] = [column]
        else
            push!(variable_names[var_name],column)
        end
    end

    variables = Dict{String, Any}()
    for (i,var_name) in enumerate(keys(variable_names))
        if length(variable_names[var_name])==1
            variables[var_name] = vec(Array(df[:,variable_names[var_name]]))
        else
            variables[var_name] = Matrix(Array(df[:,variable_names[var_name]]))
        end
    end
    return variables
end

# print("Input for Matlab:", '\n')
# print("previn = ", previn, ";", '\n')
# print("pfmin = ", pfmin, ";", '\n')
# print("pfminneg = ", pfminneg, ";", '\n')
# print("actualdiseases = ", actualdiseases, ";", '\n')
# print('\n')

# previn = [0.0230534577526959, 0.0295250422798251, 0.0247689382994323, 0.0315067109244799, 0.0370991081742557, 0.0366502056733995, 0.0279097952776079, 0.0389762037833059, 0.0208943570101915, 0.0222453349734917, 0.0283895950994135, 0.0281041742634757, 0.0342985264830265, 0.0242286689245111, 0.0300614228338476, 0.0324587582674629, 0.0331451364472201, 0.0214553620186261, 0.0376763711574193, 0.0333655926579872, 0.0406164784647051, 0.031425905301204, 0.0321342516504009, 0.0355784380213378, 0.0277670170266779, 0.0358148884628447, 0.0308255443922441, 0.0275994016686287, 0.0218527786117252, 0.0239788881066742, 0.0261490891718568, 0.0230553022691171, 0.0239360301905651, 0.0261355970001302, 0.0373705016100649, 0.0363368938436778, 0.0263577018130165, 0.0227359522657424, 0.0242710354428358, 0.0271638391160356, 0.0226312710695577, 0.0386955441559298, 0.0238067460328314, 0.0252640253833867, 0.0232042145005477, 0.0255529932593705, 0.0354832878574329, 0.040694790120068, 0.0387358719899655, 0.031159561584258, 0.0283243162115461, 0.0385320733245259, 0.0321850399295379, 0.0228186443222404, 0.0335460345184117, 0.038893794180543, 0.0268249508247191, 0.0334115239425127, 0.0303120045577703, 0.0348291910555565, 0.0294613124651552, 0.0362975070745062, 0.0246751020754088, 0.0305195759815281, 0.0272960563709025, 0.0334652915073471, 0.0259994958611438, 0.0267098335990337, 0.040647585382083, 0.0216562044170748, 0.0307418446139463, 0.0319384101710453, 0.0275240786731279, 0.0394149146520865, 0.0228506396512835, 0.0367640885629434, 0.0264953246489228, 0.037440034209939, 0.0297407339116828, 0.0385719837844939, 0.0271356024747371, 0.0285817282509067, 0.0240694693066806, 0.0358802060887939, 0.0245440462112277, 0.0305924551771928, 0.0387961661068613, 0.0235170936321763, 0.0387107007254399, 0.0213565524874163, 0.0319289788795099, 0.0251953222504976, 0.0352329663284143, 0.0230140103828934, 0.0339764229446642, 0.0288298728574788, 0.0223022346839054, 0.0274592547394084, 0.0327885358056755, 0.0260541875005938];
# pfmin = [
#     0.862391858082073 0.614624385275664 0.539611385608728 0.632831241994659 1.0 1.0 0.774720006170335 1.0 1.0 0.56244856831602 0.736447721537061 0.640988266310589 0.726667100665582 0.756995896214162 1.0 1.0 0.717634674942516 1.0 0.783285513928666 1.0 0.569119379571531 0.669258072685321 0.625532570105459 1.0 1.0 1.0 0.815644367088964 0.538314966898565 0.526408626046779 1.0 0.77442935348819 1.0 1.0 1.0 1.0 1.0 0.66918549981669 1.0 0.742330874960291 1.0 1.0 0.641775381751934 1.0 0.801233202937997 1.0 0.677794995607488 0.51091163543719 0.792015335139736 0.650494614315876 0.558506245213354 1.0 0.717597821828022 1.0 0.62060884614105 0.651086147138882 0.530511829865759 0.589212423747889 1.0 0.864956837654359 1.0 0.672696932067646 0.605238200527946 1.0 1.0 1.0 1.0 0.504092936324235 1.0 1.0 1.0 0.516778764354619 1.0 0.735984958496777 0.817961316827666 0.761610934236559 0.701282866987349 1.0 1.0 1.0 0.883971199927949 0.600212319838954 0.818222112432571 0.516038418871991 0.703392207736629 1.0 1.0 1.0 1.0 1.0 0.77771500078052 1.0 0.82171040008863 0.790495424108934 0.802622563007498 1.0 0.605058914958124 1.0 1.0 1.0 1.0
#     1.0 0.860786817985174 0.830713820997187 0.833299333016744 1.0 0.691005159495752 1.0 1.0 1.0 0.887298212389558 0.836311276586702 0.5715842827979 1.0 0.533611824530981 0.664732226334532 1.0 0.847180787623364 0.611990639562347 0.730953248028717 1.0 0.876439973812396 0.89282201719843 1.0 0.832545258787233 1.0 0.510837863005187 0.747739152074005 1.0 1.0 0.692443355290812 1.0 1.0 1.0 1.0 0.501166617700821 0.551443578293769 0.564259631683136 0.683481055252833 1.0 0.629254268444593 1.0 1.0 0.769143948076114 1.0 0.629765427713757 0.896810648246045 0.882914828811803 0.874686559225294 1.0 0.517542477109637 0.882121878562669 0.794284132020501 1.0 0.787517100794088 1.0 0.752546038999182 1.0 1.0 1.0 1.0 1.0 0.630874816487961 1.0 0.673246051100155 1.0 0.620178097070097 1.0 0.758343464326249 0.757548136411276 0.754531473667267 0.64334502493445 0.567820661805211 0.650364309704737 0.777639446719762 0.886356931288001 1.0 0.628582606357248 1.0 1.0 1.0 0.6332778332086 0.588474125299077 1.0 0.595159226778034 0.620369641102784 0.830564759595417 0.553064970080255 1.0 0.861513180073713 1.0 0.685527271176655 0.816087868263427 0.764379134816018 1.0 1.0 0.521235673019451 1.0 0.726414162318273 1.0 0.517755836536797
#     1.0 0.52355117309113 0.88051868647936 0.801647829912284 0.873371911658049 1.0 0.659388444272828 1.0 0.640142141990192 1.0 0.774999076403831 0.535777633557932 1.0 0.730477553866881 1.0 0.553403400675336 1.0 1.0 0.867932634049189 1.0 0.515713954164648 0.823145668714236 0.715919063271888 1.0 1.0 1.0 1.0 0.82090122526036 1.0 0.522219906061201 0.891535722615343 0.810236736684924 0.535930821168078 0.534500138628072 0.680370489930259 0.537081480704547 1.0 1.0 0.776697720574664 1.0 0.512347469839164 1.0 1.0 1.0 0.560943990481246 0.613204333156257 1.0 0.843567520157576 1.0 0.838697344194098 1.0 0.742290347955886 0.894006847466984 0.822739885664364 1.0 1.0 0.792612869845648 1.0 0.787843960804582 1.0 0.885162876188285 0.892040318333855 0.882287232808117 0.793330793336859 0.798304473865477 0.665886719363891 1.0 0.880278858481583 1.0 1.0 1.0 0.570241185488709 1.0 1.0 1.0 0.616907334412528 0.527239932152065 0.758121891533862 0.70563090159334 0.53435649341315 1.0 1.0 0.577936389949986 1.0 0.611015980172476 1.0 0.604739315906822 1.0 1.0 0.832158831111182 1.0 0.812653724982513 1.0 1.0 1.0 1.0 1.0 0.648949282177767 0.636973827208341 0.868780884772503
#     1.0 1.0 1.0 1.0 0.555525254996741 1.0 1.0 1.0 1.0 1.0 0.522261579052455 1.0 1.0 0.7773388527326 1.0 0.840215521773613 1.0 0.872777672234718 0.791362266129775 1.0 0.553762194412106 1.0 0.606568850442215 1.0 0.529699747758896 1.0 0.536475914550133 0.555548422555472 1.0 1.0 0.804876343921933 1.0 0.888332325254809 1.0 0.895816875766749 0.538774914756015 1.0 0.897521113909547 1.0 1.0 1.0 1.0 0.745404247464038 1.0 0.561716006705537 1.0 1.0 0.890969606917799 0.598851285929725 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.705739191606213 1.0 0.568278589630638 1.0 0.899633494660176 1.0 1.0 0.68975559086303 1.0 1.0 0.587329724805593 1.0 0.558110856370695 0.596009310919193 0.51066281322375 1.0 0.506965980403473 1.0 1.0 0.549508252033161 1.0 1.0 0.873886110920208 0.866941447531518 0.61970361689821 1.0 0.725008135485891 0.603111005176033 1.0 1.0 0.692478668601429 1.0 0.71654398342641 0.76000200079136 0.715280188009892 0.757446518463546 0.715730909534603 0.767427352677797 1.0 1.0 1.0 0.651157773273781
#     1.0 0.845179586909148 0.678800147279475 0.666947828988192 0.845134222410332 1.0 0.881258106504844 1.0 0.861909954302952 1.0 0.576073529767724 1.0 0.60044570885444 0.798944277149329 0.875146615643235 1.0 0.879898527052726 1.0 1.0 1.0 1.0 0.676171349448536 0.784805964876715 1.0 1.0 1.0 1.0 0.54279981407957 1.0 0.683837380253475 0.634593810407812 0.826464567721609 0.7642915943269 1.0 1.0 1.0 1.0 0.561070154373759 1.0 0.71667779929923 0.584623189507193 1.0 1.0 1.0 0.889575415681304 0.761293469192013 1.0 0.502824005793984 1.0 0.88013031202506 1.0 1.0 0.826365522428789 0.859902610816706 1.0 1.0 1.0 1.0 1.0 1.0 0.652065784035735 0.899292252288241 1.0 1.0 0.726720536610384 0.707346099954725 0.691859847344518 1.0 1.0 1.0 0.71276129812823 1.0 0.833805613804435 0.605076581182042 0.77941512945059 0.667503432882702 1.0 1.0 0.86378770369161 0.588638283224652 0.523686452320384 0.818209549525587 0.707739335402653 0.789806588627613 1.0 1.0 0.703127500644452 1.0 0.655879283092838 0.611079819689346 0.797064828928721 1.0 0.738647458466172 1.0 0.74549275628993 0.777189347098106 0.896076719447527 0.886697166669508 0.591892932590046 1.0
#     1.0 0.745528974399747 1.0 0.871510992747256 1.0 0.769621116233965 0.587747925721115 0.775040866216699 1.0 0.668573633442596 1.0 0.59462765695519 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.720181420401439 1.0 0.772486517769296 1.0 1.0 0.594446621999595 1.0 1.0 0.832689945876027 0.718651003547733 1.0 0.687247897959049 1.0 0.793235207515593 1.0 0.650217250557282 0.513095379167584 0.542584674045156 0.850828447782113 1.0 0.789877248786552 0.878697214639874 0.638725608090764 0.773912008957694 0.668864672699237 1.0 1.0 0.869900455904204 0.514718998713404 0.84296417853734 1.0 0.590787147416059 0.788653605888949 1.0 0.686060044546084 0.617556804312788 1.0 1.0 0.802745092497431 1.0 1.0 1.0 1.0 1.0 1.0 0.715026185016993 0.617900860943861 1.0 1.0 0.838518579322061 0.746522668208558 0.621616944300797 0.780679666059636 0.881135385288465 1.0 0.603784710679393 0.535747652449323 1.0 0.666950019212361 0.664579834348701 1.0 0.89544828723353 1.0 0.796865539863014 0.511726083890146 1.0 1.0 1.0 0.827477666982353 0.832014712331768 0.631098381249039 0.575068784897995 0.658348721324493 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
#     0.775859418873685 0.826576543772075 1.0 1.0 1.0 0.797664946003434 0.621410563425435 1.0 1.0 0.633194795310146 1.0 0.756375821439248 0.708642557742243 1.0 0.582273508058735 1.0 0.784112734935844 1.0 1.0 1.0 0.745863559924422 1.0 0.701862105036272 1.0 1.0 1.0 0.598859815002201 0.609286661850838 0.683945756557354 0.57699894133576 0.562105939180104 1.0 1.0 0.700542741019971 0.575780355017871 0.614642481324429 1.0 1.0 1.0 1.0 0.875488022820075 1.0 0.515841000636059 1.0 1.0 1.0 0.614822426719444 1.0 0.786031760263154 1.0 0.684553292568229 0.680775169447095 1.0 1.0 1.0 1.0 0.777113783158879 0.748195196245313 0.748794453754915 0.733905851363729 1.0 1.0 1.0 0.680080164033612 1.0 1.0 0.614112986868214 1.0 0.683064853984617 1.0 1.0 1.0 1.0 0.651885135683459 0.896809158550798 0.738776959730088 1.0 0.866634365688325 0.818254947257015 0.508505301516152 0.647547817300927 0.706623375704015 1.0 0.510698534182614 1.0 0.588987669935141 0.782599760839314 1.0 1.0 1.0 0.675149890014186 1.0 0.505203564652476 0.737562022933052 0.648332007439537 1.0 1.0 1.0 1.0 0.691004603504736
#     0.679415373243575 0.615708298691753 1.0 0.869275528254612 1.0 0.675680684367202 0.503704127526957 1.0 0.877303810003354 0.865863699049843 0.731937190653154 0.716940643582622 1.0 1.0 1.0 1.0 1.0 0.555494478135215 0.841035674749405 0.746853297765801 0.611611791399126 1.0 1.0 0.631837424980928 1.0 0.739103398651581 0.746147063239919 1.0 1.0 1.0 0.592592158655986 0.841659105859701 1.0 0.570818790511959 1.0 1.0 0.713458498154279 1.0 0.645631363676395 0.640919029795105 0.837557442489292 1.0 1.0 0.553932548899079 1.0 1.0 0.710051488128538 0.556593574092057 0.792557912955829 0.567197897102117 0.566060635004416 1.0 0.62531070517826 0.884197905534879 0.509538813825981 1.0 0.597015154478678 0.731559834113342 1.0 0.847918879548944 0.687612806463705 0.531602729809758 1.0 0.853258595894007 0.520723589426472 0.626879644876862 0.763615466248782 1.0 1.0 1.0 0.859540000825689 1.0 0.621331744036559 1.0 1.0 0.709627783961356 0.523301755612788 0.858997782248472 1.0 1.0 0.812426226932356 0.691904299819083 0.892004851708365 0.862756251865991 1.0 1.0 1.0 1.0 0.814130321206144 1.0 0.65036506817971 1.0 1.0 1.0 1.0 1.0 0.785387485160023 1.0 0.633839522308075 1.0
#     0.821590616450953 0.812730511801121 1.0 0.725092880737099 1.0 0.541898626507678 0.528584253719656 0.746663634751677 1.0 0.880508913128054 1.0 0.684783577111367 1.0 0.637460589923562 0.748942161388156 1.0 0.727157242350375 1.0 0.793913598312554 0.707373022795305 0.679323568446823 1.0 1.0 1.0 0.556407085861237 1.0 1.0 0.790445414482601 0.630016281525929 1.0 0.812917919914669 1.0 0.856103208065662 1.0 1.0 1.0 0.520366611973991 1.0 0.799738345736247 0.71077528865137 0.743030723618569 1.0 1.0 1.0 1.0 0.529631167380151 0.768184125863819 1.0 1.0 0.839453947046616 0.54196890851229 1.0 0.833670895386048 0.770195736490985 0.714574902873354 0.693525631645462 1.0 1.0 1.0 0.690066089198524 0.859727758271647 1.0 1.0 0.543324462270212 1.0 0.684632939171217 1.0 0.839798458651428 0.742019913386703 0.694227745884942 1.0 0.705590095388846 1.0 0.682000044408579 1.0 0.873082993575727 0.528663528117989 0.697588156044328 1.0 0.804131224492113 1.0 1.0 1.0 0.760392101111434 0.540239198610207 1.0 1.0 0.72903330478577 0.802184714446892 0.72524247496084 1.0 1.0 0.563511783443164 0.686496643620181 0.880287435914184 0.502649896011392 1.0 0.630342159905469 0.727917002730791 1.0
#     1.0 0.887738938927808 1.0 0.693351724102034 1.0 1.0 0.802565335102806 1.0 0.572502158427843 1.0 1.0 1.0 1.0 1.0 1.0 0.812958479592228 1.0 1.0 1.0 1.0 0.801119050149727 1.0 1.0 0.813197799232727 1.0 0.873925810760371 1.0 1.0 0.59312757677032 0.693594271575613 1.0 0.893176052505279 1.0 1.0 1.0 0.688461599755851 0.675950764437403 1.0 1.0 1.0 0.795093024323687 1.0 0.529473249682468 1.0 0.613319092648085 1.0 1.0 0.791333824694362 0.896368011275628 1.0 1.0 0.738166957295671 0.79889860202692 0.87552967481354 1.0 0.601733985541742 0.692626125747447 0.71517837396601 1.0 0.672142983643718 0.836453856026667 1.0 0.631775542347357 1.0 1.0 0.52390887085991 0.657885474539834 1.0 1.0 0.843943189224939 0.740299545371978 1.0 0.561309684062669 1.0 1.0 1.0 1.0 0.604304577110903 1.0 0.703540233437683 1.0 1.0 0.606791079605559 0.72502452663972 0.855321288184733 1.0 1.0 0.789298692951809 1.0 0.514526611023174 0.783172680720558 1.0 1.0 0.596955801752549 0.831071859395255 0.567199048881033 0.513403690775022 1.0 1.0 1.0
#     0.714162558489348 0.812048003686532 1.0 0.836797817393158 1.0 1.0 1.0 1.0 1.0 0.632767960297795 0.613636920868077 1.0 1.0 0.718528467989333 1.0 1.0 0.608603820164101 1.0 0.7632788702474 0.621111333426235 0.512636177778856 0.865300484133782 0.647892264714754 0.859450545978349 0.782821431306473 1.0 1.0 1.0 1.0 1.0 0.758232206995281 1.0 1.0 1.0 1.0 0.535746181901332 0.60493075046971 0.659540855614731 1.0 1.0 0.775536590958589 1.0 1.0 0.718706086337803 0.619631480145613 0.778566553374183 1.0 1.0 1.0 0.877420588389829 0.511349485976333 0.703334742132116 0.666421994792658 1.0 0.876019995403669 0.622654226436757 1.0 1.0 1.0 1.0 0.545250278036024 1.0 0.834514780415655 0.581721348039424 1.0 0.694327807489012 1.0 0.596711861085944 1.0 0.872370111080426 0.688108874554522 1.0 1.0 1.0 0.51962776450243 1.0 1.0 1.0 0.572315914904144 0.574751248597942 1.0 0.762527898723983 0.761431132203176 1.0 1.0 1.0 1.0 0.539348074936494 1.0 1.0 1.0 0.881853263406422 1.0 0.676217903927729 1.0 1.0 1.0 1.0 1.0 1.0
#     1.0 0.607411133805864 0.642348744236827 0.840057831626649 1.0 1.0 1.0 0.671122977380447 1.0 1.0 1.0 0.783747977955432 0.528027495777994 1.0 0.881467077557915 0.675740788779614 1.0 1.0 1.0 0.71458120219997 0.548068324492733 1.0 1.0 0.829851999222396 1.0 0.89009626080556 0.539795041750562 0.83663620943275 0.718910136536295 1.0 0.6954338320732 0.865850077609041 1.0 1.0 1.0 0.593941512265415 1.0 1.0 0.753910174764683 1.0 0.702083366082275 1.0 1.0 1.0 1.0 1.0 1.0 0.87295475817441 0.569971237175763 1.0 1.0 0.7537506166062 1.0 0.729048209807477 1.0 0.513246360221801 1.0 1.0 0.604254958040275 0.779321751635377 0.604788194272456 0.540415412732681 0.864059381752319 1.0 1.0 1.0 0.531815422772642 1.0 0.751439138259616 1.0 0.662409288795912 0.637985279526625 0.504308550634879 0.555483268320448 1.0 1.0 1.0 0.890646282111115 0.861002581935186 1.0 0.69163146184836 1.0 0.65082623221787 0.866305396390449 1.0 0.591510509481204 1.0 1.0 0.862928535079577 1.0 0.567879223250033 1.0 0.785978974813633 0.849552231153498 0.676518372921096 0.695281349139541 0.714384279048261 1.0 1.0 0.690504218748814
#     0.610828805241403 1.0 0.876980680707237 1.0 1.0 1.0 0.503541291519306 0.749446818931176 0.889924584336572 1.0 1.0 0.806620002751313 1.0 0.777520911034116 1.0 1.0 1.0 0.557150382698488 1.0 0.582239305404797 0.768192662826308 0.785449565649651 0.836995888339401 1.0 1.0 0.716460792160418 0.881371693427075 0.758259107499707 0.565537750649678 0.755139412838717 1.0 1.0 0.592763481556477 1.0 1.0 0.797855102337294 1.0 1.0 0.781793864195085 1.0 0.536577647253906 1.0 0.598534541440897 0.702862962763485 0.536818048683357 1.0 0.557439899497972 1.0 1.0 1.0 0.87257567268919 0.647244579459842 0.811570052457665 0.611705580974673 1.0 0.85954996192427 1.0 1.0 0.641909530244065 1.0 1.0 0.880325841663788 1.0 0.642681025432666 0.770089250181547 0.603556472198922 0.775223090525384 1.0 0.866008575197218 0.743784367368329 1.0 0.852262991554147 1.0 0.802420039335204 1.0 0.773099966489825 1.0 0.802315759708556 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.525694563327408 1.0 0.824879666019106 1.0 0.834656270900896 1.0 1.0 1.0 1.0 1.0 1.0 0.882406039060536 0.658947416181967 1.0 1.0
#     0.713129763157058 1.0 0.715943154933403 1.0 1.0 1.0 0.71676748737313 0.846614078607634 1.0 1.0 0.846412915630927 0.767863751222026 0.849555311806608 1.0 1.0 1.0 1.0 1.0 0.783335205757558 1.0 1.0 0.78599803810389 0.77099155118853 1.0 1.0 1.0 0.704718770676202 0.517129907092343 0.864469063410502 0.891337658155561 0.845422368839252 1.0 0.781622476059683 1.0 1.0 0.612954269138826 0.850488620454726 0.643295790895639 0.618844617041912 0.552252520550654 1.0 1.0 0.680342731500064 1.0 1.0 0.664289904131659 1.0 1.0 0.605525715737612 0.785459827710875 0.858171743931137 0.53675298782118 0.602046765745844 0.618196442529072 1.0 0.826505735118729 1.0 1.0 1.0 0.510058454979366 0.573244298684134 1.0 1.0 1.0 0.816584461733396 1.0 0.810831308587531 1.0 1.0 0.531420897622577 0.851449837367238 0.808064786277026 0.764345815480992 1.0 1.0 0.889346040972149 1.0 1.0 1.0 0.571262417127173 0.677936667084438 1.0 0.520458992915428 0.842671399029026 1.0 0.708844626111373 1.0 0.752955112332148 1.0 1.0 1.0 0.55979648021618 0.79167930227455 0.70646100009273 1.0 0.837157794811666 0.655009325547472 0.850558803767623 0.843241452235305 1.0
#     0.629865097185991 0.851852474500192 0.677009184577276 1.0 0.514517367260893 0.558396775651665 1.0 0.532204302179506 1.0 1.0 0.88150797119269 0.78858143004308 1.0 1.0 0.641881483317918 1.0 1.0 0.828436793204586 1.0 0.636342602933176 0.597483066239744 1.0 0.839059237325111 0.777051387967866 1.0 0.803969687625906 1.0 1.0 1.0 0.799203438420723 1.0 1.0 0.748232223236317 1.0 1.0 0.575931785392725 0.723541253164483 1.0 0.579293869575017 1.0 1.0 1.0 1.0 0.741825670104616 1.0 1.0 1.0 1.0 0.801896165543568 0.863375822976042 0.639798341222499 1.0 0.833953237166011 1.0 0.651690459499945 1.0 1.0 1.0 1.0 1.0 0.81233545658834 0.708787945461449 1.0 0.704573167455062 1.0 1.0 1.0 0.549677342770199 1.0 0.592606335872201 1.0 0.604441491583518 0.576894720359455 1.0 0.508141000028147 0.666681499793433 1.0 1.0 0.613662040766524 1.0 0.749256612474568 1.0 0.706153791600774 0.880668573591029 1.0 0.507036270767329 1.0 1.0 0.766487053573091 1.0 0.688552886951121 1.0 0.510886867992261 0.877367846254966 1.0 1.0 0.504860742253793 1.0 0.5226390347757 0.52732914942393
#     0.896634009897625 1.0 1.0 1.0 1.0 0.647964495310907 0.56331893423227 1.0 0.655617455867081 1.0 1.0 0.719477187095289 1.0 1.0 0.823003790466973 1.0 1.0 1.0 0.54417293550497 0.596069886193695 0.844211000153263 1.0 1.0 0.529888060176286 1.0 0.708609056486596 0.884229668837625 0.790949155510243 0.852276943544415 0.76425074861096 1.0 1.0 1.0 1.0 0.721663733314736 1.0 1.0 1.0 0.665716236079822 1.0 1.0 0.713200944366519 0.615091548171175 0.699044796620894 1.0 1.0 1.0 1.0 0.713701927838012 1.0 1.0 1.0 1.0 1.0 0.634102103474393 1.0 1.0 1.0 0.664573484661703 1.0 1.0 0.821993850192908 0.847412255511799 1.0 1.0 0.563870994267807 0.709639841226599 0.844028638327023 1.0 1.0 0.636193953991219 1.0 1.0 1.0 0.698952740792528 0.535529064163145 1.0 0.652005087942678 1.0 0.512223414838388 0.795408580179716 0.738893351820224 1.0 1.0 0.500283707576383 1.0 1.0 1.0 1.0 0.667652593466797 1.0 1.0 0.863882187496974 0.68176708431576 0.58908315767052 1.0 1.0 1.0 1.0 0.815037070240068
#     0.811590622551033 1.0 1.0 0.630477160617141 1.0 1.0 0.560909656773907 1.0 0.747283588936391 1.0 1.0 0.689390974245193 1.0 1.0 1.0 1.0 0.898000569268923 0.762278197579971 1.0 1.0 0.652623875310333 1.0 0.602343885565359 0.787865560939045 1.0 0.522993444641918 0.854177719615813 1.0 0.733304838405505 1.0 0.507890071594357 0.753570462304003 0.64481499045996 0.825239730082439 0.557086084567814 0.525175914242809 1.0 0.889638548169369 1.0 1.0 1.0 0.519177852135886 1.0 0.601386358094721 1.0 1.0 1.0 0.613839102344973 0.848599414323864 1.0 0.87956786502243 0.511747569277958 1.0 1.0 0.67195172859619 0.632555288116813 0.567461262578798 1.0 1.0 0.640173477322355 1.0 0.513936411275667 0.854612282843894 0.551790512025495 0.786243121922089 0.759281192156933 0.828045067291944 1.0 0.705069152561263 1.0 0.679128443722517 0.772713182815157 1.0 0.794259996012433 1.0 1.0 1.0 0.633574118511481 1.0 1.0 0.838926532771428 1.0 1.0 0.827646911555992 1.0 1.0 1.0 1.0 1.0 1.0 0.742024069399952 0.810654055418215 1.0 1.0 1.0 0.766330308330617 1.0 1.0 1.0 0.844829421294637
#     0.801109359807527 1.0 1.0 1.0 0.627352509254136 1.0 0.835434661548468 0.705835419624546 0.582273394969775 0.719429715956735 1.0 1.0 0.643481767124917 0.58261377365715 1.0 1.0 1.0 1.0 1.0 1.0 0.793667047775803 1.0 1.0 1.0 0.73706619588462 0.861664446588489 0.791314372583316 0.694442690671651 0.623866504652936 0.802501397157775 1.0 0.513264250777693 0.705134642400637 1.0 0.781535272735231 1.0 1.0 0.505594210201475 0.820615530346344 0.669996906896014 0.849998019013564 1.0 0.587628146919038 1.0 0.572139501221982 0.710417106540087 0.717316477282088 0.573054234414607 1.0 1.0 1.0 0.819867879259275 0.806557630301688 0.875307521329996 1.0 1.0 1.0 0.611461161928043 1.0 1.0 1.0 0.577961498214857 1.0 1.0 1.0 1.0 1.0 0.833775178324811 1.0 0.827986640882024 1.0 0.815882737808481 1.0 1.0 0.79800919290587 1.0 0.81058674314451 1.0 0.80178694185544 0.592212478528695 0.781767279247413 1.0 0.607013240408387 1.0 0.643335752147823 0.500206233626354 0.675759825992872 1.0 1.0 0.724411653620926 0.635184278201717 1.0 1.0 1.0 1.0 0.829075205627433 0.629330973215722 1.0 0.548982080391441 0.69064332279014
#     0.529948232719195 0.749836836390569 0.846605888264962 1.0 1.0 0.795141654885818 0.813112243455962 1.0 0.837986826281153 1.0 1.0 1.0 1.0 1.0 0.661243702930808 0.810025007128282 0.891065406043989 1.0 0.700254541307114 0.636038384927918 0.894819635730166 0.679132813856353 1.0 0.773647064105544 1.0 1.0 1.0 1.0 0.517966901041099 0.606039385726497 0.722336606078171 0.610661626297958 1.0 0.586077516070982 1.0 0.596681935339553 1.0 1.0 1.0 1.0 0.604359266237349 0.879269491708588 1.0 0.874496043492226 1.0 0.620499311092207 1.0 0.658176022297171 0.526712546222357 0.570402849714304 0.722330207196057 1.0 1.0 0.59488017223918 1.0 1.0 0.574376899728744 0.565488826927913 0.527962498577524 1.0 1.0 1.0 0.593843716433811 1.0 0.59063926061402 0.835147687741612 1.0 0.547031867537728 0.529457916560545 1.0 0.797982478116287 1.0 1.0 1.0 0.87146052696372 1.0 0.853408501915253 1.0 1.0 0.732599539357147 0.671463779935659 1.0 0.747084862055371 1.0 0.584859875038807 0.598086364510292 1.0 0.798003022118585 0.65793605775236 1.0 0.694468835055477 1.0 1.0 0.847645621493114 1.0 0.628123846811497 0.883417326043936 0.506571637771072 1.0 1.0
#     1.0 1.0 1.0 0.665231844478601 0.631872749351207 1.0 0.531512866221767 0.76669011508669 0.50239537602979 0.725146870661367 0.843487786737933 0.821376994488768 0.714957106316158 1.0 1.0 0.507917316784404 0.564767387389055 0.721332941269063 0.727505988138869 0.828017884035321 1.0 0.74076862755572 1.0 1.0 0.795194043171647 0.753002255942757 1.0 1.0 1.0 0.706145454535693 0.688345809701489 0.596055668110888 0.850668895713487 0.732294144765905 1.0 0.771056774442598 1.0 0.778158677101542 1.0 0.74129252697546 0.537603100831477 1.0 1.0 0.568567968500155 1.0 0.583263282606263 1.0 1.0 0.731000255220414 0.66049844911258 0.721702758704671 1.0 0.899562796513289 0.519397601124674 0.74394406594035 1.0 1.0 0.70863313699177 0.598902096800638 0.53772452425603 0.838084141047979 0.873556772578828 1.0 0.759134010562808 1.0 1.0 1.0 1.0 0.721783619465259 0.512318841803552 1.0 0.851582473408674 1.0 0.526699719946876 1.0 0.698075206847683 0.624902620269432 1.0 0.656332947522915 1.0 1.0 1.0 1.0 1.0 0.827453103681166 0.82480238189912 0.835617925939688 1.0 1.0 1.0 0.708936802079286 0.514309701323157 0.566592908013167 1.0 1.0 0.719819033365586 0.887694729339644 0.665150623621201 1.0 0.559252105198184
#     1.0 0.719994907222631 0.890389418006404 0.784718074413879 1.0 1.0 1.0 0.785386163433361 0.692869416083972 1.0 0.835637977827095 1.0 1.0 1.0 0.558595906601859 0.574144281428735 0.757727029541309 0.708602062206539 1.0 1.0 0.797177601190916 1.0 0.586716604870966 1.0 1.0 1.0 0.516803184337581 0.753841019453524 1.0 0.818808995112029 0.882314581666575 0.883434615115523 1.0 1.0 0.568151000038167 0.728792422964076 1.0 1.0 0.823103082699972 0.670287536629337 1.0 0.691066326500884 0.606421036229025 1.0 0.730477234424218 1.0 0.528143538412895 0.613696652673074 1.0 1.0 0.725272684037525 0.811620160202366 1.0 0.888361315135842 0.756966452731424 0.845304858813831 1.0 1.0 0.831361477888759 0.887199591671565 0.860759552166712 1.0 1.0 1.0 0.613137278667072 1.0 1.0 0.795017251806787 0.719547318476881 0.880793558458131 1.0 0.674330283134515 0.726332470787075 0.586053802470759 0.86057933471753 1.0 1.0 0.519020323768843 1.0 1.0 1.0 0.528190931592716 1.0 0.833900722112758 1.0 1.0 0.565153400101294 1.0 1.0 0.612361941437756 0.594925860647369 1.0 1.0 1.0 0.68224762651082 1.0 1.0 0.735135151742215 0.741088182344272 1.0
#     0.704764137910723 0.516543060009118 0.787363129010562 1.0 1.0 1.0 0.787684816330388 0.809021742626919 0.858940492381161 0.813314856198714 0.854982219698155 0.565464008863519 0.603308871199958 1.0 1.0 1.0 0.743636875019999 1.0 0.690712083564766 1.0 0.538799448469953 0.861828362354903 0.856293844864486 0.839042691126115 0.617339965467687 1.0 1.0 1.0 0.805825201922403 1.0 1.0 0.652656333327164 1.0 1.0 1.0 1.0 0.791770566587427 0.733190111399159 0.824429129711086 0.822382523829836 0.634954360431385 1.0 1.0 1.0 0.600379980365539 1.0 0.61341017675448 0.734032081055431 1.0 1.0 1.0 0.648083384326814 1.0 0.500417250725503 0.570776552679572 1.0 0.75304604527309 1.0 0.771754580840548 1.0 0.53511758916726 0.538871747816814 1.0 1.0 1.0 1.0 1.0 0.863729891809509 1.0 1.0 0.749101646676876 1.0 0.613349108245513 0.639945619414975 1.0 1.0 0.87045733331141 1.0 0.81271336278424 0.731509167323403 0.612679926353741 0.798316874291633 0.833496661619177 0.534662551114966 0.529844143278882 1.0 0.806392611861079 0.892357142035264 0.650396122880931 1.0 0.846986338769314 1.0 1.0 1.0 0.854439691498227 0.894988420682706 0.899611137021419 1.0 0.585412012488945 0.559893812513221
# ]
# pfminneg = [1.54267251412907e-6, 1.17964085860106e-6, 2.68412323751408e-7, 5.42528707960789e-7, 4.5919663822241e-6, 5.08433659894025e-8, 4.17476330040942e-6, 4.64062877198977e-8, 1.14162664641309e-6, 1.11650783823763e-7, 3.05749294614222e-7, 3.93509997309251e-7, 4.51127927072701e-9, 7.75436961895553e-7, 9.7403117595162e-8, 8.79497010620158e-8, 1.54742494945183e-7, 1.90579886024137e-9, 2.71922286780485e-6, 2.54475271858773e-7, 1.10508971669483e-6, 3.5275867403085e-8, 2.7888664751707e-7, 2.14233736269678e-6, 4.26554155750906e-7, 1.02932460121662e-6, 3.41943729551841e-7, 1.29605903707596e-6, 1.33838096056008e-6, 3.46641697292869e-7, 1.73613421871378e-6, 3.29600615285585e-6, 4.04857320365692e-7, 2.18667236399945e-7, 4.07746597808639e-7, 2.62734044449243e-7, 3.78236019179211e-7, 8.29908724141745e-6, 8.41841851033033e-8, 2.36056877111509e-5, 4.8968307584697e-6, 1.18321836741017e-8, 1.14760718369408e-7, 1.71921006662562e-7, 3.02025358553053e-6, 3.28081338165406e-7, 6.6009297290898e-9, 1.39748716764748e-6, 4.69408065275208e-6, 1.48889608852282e-6, 1.28985905909327e-5, 1.22395603014256e-6, 7.09096693033978e-7, 1.05848765485788e-5, 4.90704955728268e-7, 1.77010642011369e-6, 2.63579996181105e-7, 9.56455147560781e-7, 3.88235864068703e-6, 4.43655384848339e-8, 2.45675803217384e-6, 7.30320506656818e-7, 1.34584419182937e-8, 1.95373923484731e-7, 1.44649226394737e-7, 3.52678191175941e-6, 1.45836680090812e-7, 1.53859817097567e-8, 6.04449604126624e-7, 4.11418430427298e-7, 8.00652301077582e-7, 5.49301461478022e-8, 3.2333764295509e-8, 2.49843497778284e-6, 3.25209066525445e-8, 4.88540597851747e-8, 3.37138373691112e-7, 6.17007666623924e-6, 1.89448937932998e-7, 3.45422058982787e-7, 1.99904824901894e-6, 7.99755824934153e-7, 2.49079262390213e-6, 2.51933719089464e-6, 1.73011019040113e-7, 3.0080965925393e-8, 3.63517465054187e-8, 3.63839155418094e-5, 2.16860306167482e-7, 1.78430279823613e-5, 5.65410576998019e-7, 2.65714158529717e-6, 5.21265299477344e-7, 2.71928278668943e-7, 1.69615811586942e-7, 1.06504831549599e-6, 9.32968765295201e-9, 4.9828502073443e-7, 4.62486104207813e-7, 1.03065598943152e-5];
# actualdiseases = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];