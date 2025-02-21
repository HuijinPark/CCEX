#%% Print dictionary type
    ####################################################################
def print_dict(dic):
	for key in dic.keys():
		newline=False
		if isinstance(dic[key],list):
			if isinstance(dic[key][0],list):
				if len(dic[key]) > 5:
					newline=True

		if newline == True:
			for i in range(len(dic[key])):
				if i==0: print("\t",key,"(#:%d)"%len(dic[key]),"\t:")
				print("\t",dic[key][i])
		else:
			print("\t",key, "\t: ",dic[key])
