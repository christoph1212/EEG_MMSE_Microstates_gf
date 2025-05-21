function [EEG,string,coeff] = filter_silent(EEG, locutoff, hicutoff)
	[~,EEG,string,coeff] = evalc("pop_eegfiltnew(EEG,locutoff, hicutoff,[],logical(1),[],0)");
end