function input_networks=load_data(networks,method)

data = load(['./2. DREAM/DREAM Networks/network_',num2str(networks),'_',method,'.mat']);
input_networks = data.input_network;

end