function [FiberInfo] = ExtractPhInfo_All(ctrlmovlist,bapnmovlist,fakmovlist)

FiberInfo.phfak = ExtractPhInt(fakmovlist,'ph');
FiberInfo.phbapn = ExtractPhInt(bapnmovlist,'ph');
FiberInfo.phctrl = ExtractPhInt(ctrlmovlist,'ph');
FiberInfo.ectrl = ExtractPhInt(ctrlmovlist,'ecad');
FiberInfo.ebapn = ExtractPhInt(bapnmovlist,'ecad');
FiberInfo.efak = ExtractPhInt(fakmovlist,'ecad');

end