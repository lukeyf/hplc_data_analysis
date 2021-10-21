from hplc_data_anal.backend.analysis_method import *
from statistics import mean
# # # print(experimentally_monitored_data('/Users/luke/PycharmProjects/hplc_data_anal/gm_data',plot=True))
x,y= extract_data('/Users/luke/PycharmProjects/hplc_data_anal/gm_data/006-1-rxn 0C.D',wavelength_nm=310)
a,b = extract_data('/Users/luke/PycharmProjects/hplc_data_anal/gm_data/005-1-rxn 0C.D',wavelength_nm=310)
# x=[mode(y)] * len(y)
plt.plot(x,y,a,b)
plt.xlim([8,10])
plt.legend(['y','b'])
plt.show()
print(peak_properties([x,[mean(y)] * len(y)],[x,y],plot=True,plot_range=[6,9.5]))

# print(draft_reaction_score('/Users/luke/PycharmProjects/hplc_data_anal/gm_data/005-1-rxn 0C.D',[8,9],310),
#       draft_reaction_score('/Users/luke/PycharmProjects/hplc_data_anal/gm_data/006-1-rxn 0C.D',[8,9],310),
#       draft_reaction_score('/Users/luke/PycharmProjects/hplc_data_anal/gm_data/007-1-rxn 0C.D',[8,9],310),
#       draft_reaction_score('/Users/luke/PycharmProjects/hplc_data_anal/gm_data/008-1-rxn 0C.D',[8,9],310))

print(('/Users/luke/PycharmProjects/hplc_data_anal/gm_data'))