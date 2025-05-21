function [EEG, string] = resample_silent(EEG, Freq)
	[string, EEG] = evalc("pop_resample(EEG,Freq)");
end


function [EEG,string,coeff] = filter_silent(EEG, locutoff, hicutoff)
	[string, EEG] = evalc("pop_eegfiltnew(EEG,locutoff, hicutoff,[],logical(1),[],0)");
end