All = [125 123 135 156 146 161 90 190 132 161 136 165 146 148 149 149 167 186 196]; % Times of all pilots
AllOther = [135 156 146 161 132 161 136 165 146 148 149 149 167 186 196]; % Times of all minus Brain Tweakers
BestFinalA = [90 190 123 125 135 156 146 161]; % Times of finalists (final + qualifiers)
BestOfAll = [123 135 146 90 132 136 146 149 167 186 196]; % Times of the best performance of all

sort(AllOther,'ascend')
disp(['All other = ' num2str(mean(AllOther)) '+/-' num2str(std2(AllOther)) ' , N=' num2str(length(AllOther))]);
disp(['All = ' num2str(mean(All)) '+/-' num2str(std2(All)) ' , N=' num2str(length(All))]);
disp(['BestFinalA = ' num2str(mean(BestFinalA)) '+/-' num2str(std2(BestFinalA)) ' , N=' num2str(length(BestFinalA))]);
disp(['BestOfAll = ' num2str(mean(BestOfAll)) '+/-' num2str(std2(BestOfAll)) ' , N=' num2str(length(BestOfAll))]);
