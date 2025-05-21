function EEG = runica_silent(EEG, rank)
	[~, EEG] = evalc("pop_runica(EEG,'icatype','runica','pca',rank,'extended',1,'interupt','off')");
end

