import random 

def init(of,st):
	inf = open ('top'+st+'.txt','r')
	s = inf.read()
	of.write(s)
def end(of,st=''):
	inf = open ('bottom'+st+'.txt','r')
	s = inf.read()
	of.write(s)
def set_pos(of,node,node_coord):
	
	of.write('Ptr<ListPositionAllocator> pos'+str(node)+' = CreateObject<ListPositionAllocator> ();\n')
  	of.write('pos'+str(node)+'->Add (Vector ('+str(node_coord[0])+','+str(node_coord[1])+',  0.0));\n')
  	of.write('MobilityHelper mob'+str(node)+';\n')
  	of.write('mob'+str(node)+'.SetMobilityModel ("ns3::ConstantPositionMobilityModel");\n')
  	of.write('mob'+str(node)+'.SetPositionAllocator (pos'+str(node)+');\n')
  	of.write('mob'+str(node)+'.Install ('+str(node)+');\n')

def gen_enb(of,enb,enb_coord):
	# channel assignment must be done before it '+str(enb)+'
	of.write('NodeContainer enb'+str(enb)+';\n')
	of.write('enb'+str(enb)+'.Create (1);\n')
	set_pos(of,'enb'+str(enb),enb_coord);
	of.write('NetDeviceContainer enbDevs'+str(enb)+';\n')
  	of.write('enbDevs'+str(enb)+' = lteHelper->InstallEnbDevice (enb'+str(enb)+');\n')

	of.write('std::cout <<"Cell ID: " << enb'+str(enb)+'.Get (0)->GetDevice (0)->GetObject <LteEnbNetDevice> ()->GetCellId () << "(" << enb'+str(enb)+'.Get (0)->GetObject<MobilityModel> ()->GetPosition () .x << "," << enb'+str(enb)+'.Get (0)->GetObject<MobilityModel> ()->GetPosition () .y << ")" << std::endl;\n')


def gen_ue(of,ue,ue_coord,udp):
	# channel assignment must be done before it '+str(enb)+'
	of.write('NodeContainer ue'+str(ue)+';\n')
	of.write('ue'+str(ue)+'.Create (1);\n')
	set_pos(of,'ue'+str(ue),ue_coord);
	of.write('NetDeviceContainer ueDevs'+str(ue)+';\n')
  	of.write('ueDevs'+str(ue)+' = lteHelper->InstallUeDevice (ue'+str(ue)+');\n')

	if (udp):
		of.write('internet.Install (ue'+str(ue)+');\n')
	  	of.write('Ipv4InterfaceContainer ueIpIface'+str(ue)+';\n')
	  	of.write('ueIpIface'+str(ue)+' = epcHelper->AssignUeIpv4Address (NetDeviceContainer (ueDevs'+str(ue)+'));// Assign IP address to UEs, and install applications\n')


	of.write('std::cout <<"IMSI: " << ue'+str(ue)+'.Get (0)->GetDevice (0)->GetObject <LteUeNetDevice> ()->GetImsi () << "(" << ue'+str(ue)+'.Get (0)->GetObject<MobilityModel> ()->GetPosition () .x << "," << ue'+str(ue)+'.Get (0)->GetObject<MobilityModel> ()->GetPosition () .y << ")" << std::endl;\n')
	

def setActiveChannels(of,flag,Alloc):
	of.write('if (Isolated == '+str(flag)+') {\n')
	of.write('ctrl.clear();\n')
	of.write('ctrl.resize (bandwidth, false); \n')
	of.write('dl.clear(); \n')
  	of.write('dl.resize (bandwidth / rbgSize, false);  \n')


	of.write('for (int i = 0; i < bandwidth / rbgSize; i++){\n')
	of.write('if (')
	if (len(Alloc) == 0):
		of.write('false')
	for i in range(len(Alloc)):
		if (i != 0):
			of.write(' || ')
		of.write('(i >= '+str(Alloc[i][0])+'&& i <= '+str(Alloc[i][1])+')')
	of.write('){\n')
	#of.write('ctrl.at(i * rbgSize) = false;  ctrl.at(i * rbgSize + 1) = false;   // CTRL Isolation\n')

	of.write('for (int j=0; j<rbgSize; j++)\n')
	of.write('ctrl.at(i * rbgSize + j) = false; // CTRL Isolation\n')
	of.write('dl.at(i) = false;\n')
	of.write('}\n')
	of.write('else {\n')
	#of.write('ctrl.at(i * rbgSize) = true;  ctrl.at(i * rbgSize + 1) = true;   // CTRL Isolation\n')
	of.write('for (int j=0; j<rbgSize; j++)\n')
	of.write('ctrl.at(i * rbgSize + j) = true; // CTRL Isolation\n')
	of.write('dl.at(i) = true;\n')
	of.write('}\n')
	of.write('}\n')
	
	of.write('}\n')

def assign_subchannels(of,enb,Assignment_1,Assignment_2,Assignment_3,Assignment_4,share_1,share_2,share_3):
	setActiveChannels(of,1,Assignment_1)

	setActiveChannels(of,2,Assignment_2)
	
	setActiveChannels(of,3,Assignment_3)

	setActiveChannels(of,4,Assignment_4)
	
	of.write('if (Isolated == 5) {\n')
	of.write('lteHelper->SetAllocation ('+str(share_1)+');\n')
	of.write('lteHelper->setIQhopping (true);\n')
	
	of.write('}\n')

	of.write('if (Isolated == 6) {\n')
	of.write('lteHelper->SetAllocation ('+str(share_2)+');\n')
	of.write('lteHelper->setIQhopping (true);\n')
	
	of.write('}\n')

	of.write('if (Isolated == 7) {\n')
	of.write('lteHelper->SetAllocation ('+str(share_3)+');\n')
	of.write('lteHelper->setIQhopping (true);\n')
	
	of.write('}\n')

	of.write('lteHelper->setIQparams (cqi_thresh, backoff, incrementVal);\n')
	of.write('lteHelper->SetPdcchMaps (ctrl);\n')
  	of.write('lteHelper->SetDlRbgMaps (dl);\n')
  	of.write('\n')
  	of.write('\n')

def generate_topo_file(outputDir, bw,rbgsize,enbs,ues,Assignment_1,Assignment_2,Assignment_3,Assignment_4,share_1,share_2,share_3):
	outfile = open(outputDir + 'topo.cc','w')
	init(outfile,'')	
	outfile.write('uint8_t bandwidth = '+str(bw)+'; //6:15:25:50:75:100:\n')
  	outfile.write('lteHelper->SetEnbDeviceAttribute ("DlBandwidth", UintegerValue (bandwidth));\n')
  	outfile.write('lteHelper->SetEnbDeviceAttribute ("UlBandwidth", UintegerValue (bandwidth));\n')

  	outfile.write('lteHelper->SetSchedulerType ("ns3::PfFfCellFiScheduler");\n')
  	outfile.write('lteHelper->SetFfrAlgorithmType ("ns3::LteFfrCentralizedAlgorithm");\n')

  	outfile.write('uint8_t rbgSize = '+str(rbgsize)+';  // for bandwidths 1,2,2,3,4,4\n')
  	outfile.write('std::vector <bool> dl, ul, ctrl;\n')
  	outfile.write('ctrl.resize (bandwidth, false);  // false for control means, the RB is enabled\n')
  	outfile.write('dl.resize (bandwidth / rbgSize, false);  // false for dl means, the RB is enabled\n')
  	outfile.write('ul.resize (bandwidth , false); // false for ul means, the RB is enabled\n')
  	outfile.write('lteHelper->SetUlRbgMaps (ul);\n')
	outfile.write('enum EpsBearer::Qci q = EpsBearer::GBR_CONV_VOICE;\n')
  	outfile.write('EpsBearer bearer (q);\n')
	for i in range(len(enbs)):
		assign_subchannels(outfile,i,Assignment_1[i],Assignment_2[i],Assignment_3[i],Assignment_4[i],share_1[i],share_2[i],share_3[i])
		gen_enb(outfile,i,enbs[i])
		for j in range(len(ues[i])):
			gen_ue(outfile,str(i)+'_'+str(j),ues[i][j],False)
			outfile.write('lteHelper->Attach (ueDevs'+str(str(i)+'_'+str(j))+', enbDevs'+str(i)+'.Get (0));\n')
			outfile.write('lteHelper->ActivateDataRadioBearer (ueDevs'+str(str(i)+'_'+str(j))+', bearer);\n')
	end(outfile)
	outfile.close()

def makeList(time,enbs,ues):
	ret = {}
	for i in range(len(enbs)):
		ret[i] = []
		for j in range(len(ues[i])):
			ret[i].append(time)
	return ret

def generate_topo_file_udp_flows(outputDir, bw,rbgsize,enbs,ues,Assignment_1,Assignment_2,Assignment_3,Assignment_4,share_1,share_2,share_3,onTime,offTime,edges,ctrlIntfrnc):
	outfile = open(outputDir + 'topo.cc','w')

	enbDevs = {}
	ueDevs = {}

		
	stT = 0.1
	endT = 1.0
	startTime = [0.1]
	endTime = [1.1]
	onTime = [0.1]
	offTime = [0.1]

	if (len(startTime) == 1):
		startTime = makeList(startTime[0],enbs,ues)
	if (len(endTime) == 1):
		endTime = makeList(endTime[0],enbs,ues)
	if (len(onTime) == 1):
		onTime = makeList(onTime[0],enbs,ues)
	if (len(offTime) == 1):
		offTime = makeList(offTime[0],enbs,ues)
	
	init(outfile,'-udp') # static traffic generation	
	#init(outfile,'-udp-dyn') # dynamic traffic generation	
	outfile.write('uint8_t bandwidth = '+str(bw)+'; //6:15:25:50:75:100:\n')
  	outfile.write('lteHelper->SetEnbDeviceAttribute ("DlBandwidth", UintegerValue (bandwidth));\n')
  	outfile.write('lteHelper->SetEnbDeviceAttribute ("UlBandwidth", UintegerValue (bandwidth));\n')

  	outfile.write('lteHelper->SetSchedulerType ("ns3::PfFfCellFiScheduler");\n')
  	outfile.write('lteHelper->SetFfrAlgorithmType ("ns3::LteFfrCentralizedAlgorithm");\n')

  	outfile.write('InternetStackHelper internet;\n')
  	outfile.write('Ipv4StaticRoutingHelper ipv4RoutingHelper;\n')
  	outfile.write('PointToPointHelper p2ph;\n')
  	outfile.write('Ptr<Node> remoteHost = createRemoteHost(epcHelper,internet,ipv4RoutingHelper,p2ph);\n')

  	outfile.write('uint8_t rbgSize = '+str(rbgsize)+';  // for bandwidths 1,2,2,3,4,4\n')
  	outfile.write('std::vector <bool> dl, ul, ctrl;\n')
  	outfile.write('ctrl.resize (bandwidth, false);  // false for control means, the RB is enabled\n')
  	outfile.write('dl.resize (bandwidth / rbgSize, false);  // false for dl means, the RB is enabled\n')
  	outfile.write('ul.resize (bandwidth , false); // false for ul means, the RB is enabled\n')
  	outfile.write('lteHelper->SetUlRbgMaps (ul);\n')
	#outfile.write('enum EpsBearer::Qci q = EpsBearer::GBR_CONV_VOICE;\n')
  	#outfile.write('EpsBearer bearer (q);\n')
	
	for i in range(len(enbs)):
		#startTime = 0
		#inc = 0
		assign_subchannels(outfile,i,Assignment_1[i],Assignment_2[i],Assignment_3[i],Assignment_4[i],share_1[i],share_2[i],share_3[i])
		gen_enb(outfile,i,enbs[i])

		enbDevs[enbs[i]] = 'enbDevs'+str(i)

		for j in range(len(ues[i])):
			numb = random.randint(1,55)
			#outfile.write('lteHelper->SetFadingModelAttribute ("TraceFilename", StringValue ("../src/lte/model/fading-traces/traces/single'+str(numb)+'.fad"));\n')

			ueDevs[ues[i][j]] = 'ueDevs'+str(str(i)+'_'+str(j))

			gen_ue(outfile,str(i)+'_'+str(j),ues[i][j],True)
			outfile.write('lteHelper->Attach (ueDevs'+str(str(i)+'_'+str(j))+', enbDevs'+str(i)+'.Get (0));\n')
			#outfile.write('lteHelper->ActivateDataRadioBearer (ueDevs'+str(str(i)+'_'+str(j))+', bearer);\n')
			'''
			if (j < len(ues[i])/2):
				stT = 0.1
				endT = 1.0
			else:
				stT = 1.0
				endT = 2.0
			'''

			outfile.write('AppInstaller(lteHelper,ueDevs'+str(str(i)+'_'+str(j))+',enbDevs'+str(i)+',ue'+str(i)+'_'+str(j)+',epcHelper,ipv4RoutingHelper,remoteHost, 9,'+str(stT)+','+str(endT)+','+str(onTime[i][j])+','+str(offTime[i][j])+','+str(i)+','+str(j)+','+str(ctrlIntfrnc[ues[i][j]])+',70000000);\n')
			#inc += 0.1
	'''	
	#print edges
	for i in range (len(enbs)):
		node = enbs[i]
		for user in ues[i]:
			outfile.write('addMyClient('+ueDevs[user]+','+ enbDevs[node]+');\n')
		for user in edges[i]:
			outfile.write('addInterferingClient('+ueDevs[user]+','+ enbDevs[node]+');\n')
	'''
	#outfile.write('fout.close()\n;')
	end(outfile)
	outfile.close()

def generate_topo_file_wifi_flows(outputDir,enbs,ues):
	outfile = open(outputDir + 'WiFitopo.cc','w')

	rate = 70000000
	start = 0.5
	endT = 3.5
	init(outfile,'-wifi')	# static traffic generation
	#init(outfile,'-wifi-dyn')	# dynamic traffic generation
	
	for i in range(len(enbs)):
		address1 = '10.1.'+str(2*i+1)+'.0'
		address2 = '10.1.'+str(2*i+2)+'.0'
		outfile.write('NodeContainer remoteHost'+str(i)+';\n')
        	outfile.write('//remoteHost'+str(i)+'.Create (1);\n')
        	outfile.write('//stack.Install (remoteHost'+str(i)+');\n')

        	outfile.write('Ssid ssid'+str(i)+' = Ssid ("ns-3-ssid'+str(i)+'");\n')

        	outfile.write('NodeContainer wifiApNode'+str(i)+';\n')
        	outfile.write('wifiApNode'+str(i)+'.Create (1);\n')
        	outfile.write('NetDeviceContainer apDevice'+str(i)+';\n')
        	outfile.write('apDevice'+str(i)+' = setAP (wifiApNode'+str(i)+', "'+address1+'", "'+address2+'", remoteHost'+str(i)+', pointToPoint, stack, address, channel, phy, wifi, mac, ssid'+str(i)+', ipv4RoutingHelper, '+str(enbs[i][0])+', '+str(enbs[i][1])+', 1);\n')

		outfile.write('NetDeviceContainer clientDevices'+str(i)+';\n')

		for j in range(len(ues[i])):
			outfile.write('NodeContainer wifiStaNodes'+str(i)+'_'+str(j)+';\n')
        		outfile.write('wifiStaNodes'+str(i)+'_'+str(j)+'.Create (1);\n')
        
        		outfile.write('NetDeviceContainer staDevice'+str(i)+'_'+str(j)+' = setSta(wifiStaNodes'+str(i)+'_'+str(j)+', stack, channel, phy, wifi, mac, ssid'+str(i)+', '+str(ues[i][j][0])+', '+str(ues[i][j][1])+');\n')
        		outfile.write('clientDevices'+str(i)+'.Add(staDevice'+str(i)+'_'+str(j)+'.Get(0));\n')


		outfile.write('address.SetBase ("'+address2+'", "255.255.255.0");\n')
		outfile.write('Ipv4InterfaceContainer staInterface'+str(i)+' = address.Assign (clientDevices'+str(i)+');\n')
		outfile.write('Ipv4InterfaceContainer apInterface'+str(i)+' = address.Assign (apDevice'+str(i)+');\n')

		for j in range(len(ues[i])):
        		outfile.write('addWirelessRoutes(wifiStaNodes'+str(i)+'_'+str(j)+',apInterface'+str(i)+'.GetAddress (0));\n')
        		outfile.write('installApp(staInterface'+str(i)+'.GetAddress ('+str(j)+'), wifiStaNodes'+str(i)+'_'+str(j)+'.Get(0), wifiApNode'+str(i)+','+str(j)+','+str(i)+','+str(start)+','+str(endT)+','+str(rate)+');\n')

	end(outfile,'-wifi')
	outfile.close()

import math
def dist(a,b):
	x = a[0] - b[0]
	y = a[1] - b[1]
	return math.hypot(x,y)

def angle(a,b):
	x = a[0] - b[0]
	y = a[1] - b[1]
	return math.atan2(y,x)

def generate_topo_file_wifiac_flows(outputDir,enbs,ues):
	outfile = open(outputDir + 'WiFiactopo.cc','w')

	scale_factor = 50.0/290
	rate = 70000000
	start = 0.5
	endT = 3.5
	init(outfile,'-wifiac')	# static traffic generation
	#init(outfile,'-wifi-dyn')	# dynamic traffic generation


	for i in range(len(enbs)):
		for j in range(len(ues[i])):
			ang = angle(ues[i][j],enbs[i])	
			distance = dist(ues[i][j],enbs[i])	
			x = enbs[i][0] + scale_factor*distance*math.cos(ang)
			y = enbs[i][1] + scale_factor*distance*math.sin(ang)
			ues[i][j] = (x,y)

	
	for i in range(len(enbs)):
		address1 = '10.1.'+str(2*i+1)+'.0'
		address2 = '10.1.'+str(2*i+2)+'.0'
		outfile.write('NodeContainer remoteHost'+str(i)+';\n')
        	outfile.write('//remoteHost'+str(i)+'.Create (1);\n')
        	outfile.write('//stack.Install (remoteHost'+str(i)+');\n')

        	outfile.write('Ssid ssid'+str(i)+' = Ssid ("ns-3-ssid'+str(i)+'");\n')

        	outfile.write('NodeContainer wifiApNode'+str(i)+';\n')
        	outfile.write('wifiApNode'+str(i)+'.Create (1);\n')
        	outfile.write('NetDeviceContainer apDevice'+str(i)+';\n')
        	outfile.write('apDevice'+str(i)+' = setAP (wifiApNode'+str(i)+', "'+address1+'", "'+address2+'", remoteHost'+str(i)+', pointToPoint, stack, address, channel, phy, wifi, mac, ssid'+str(i)+', ipv4RoutingHelper, '+str(enbs[i][0])+', '+str(enbs[i][1])+', 1);\n')

		outfile.write('NetDeviceContainer clientDevices'+str(i)+';\n')

		for j in range(len(ues[i])):
			outfile.write('NodeContainer wifiStaNodes'+str(i)+'_'+str(j)+';\n')
        		outfile.write('wifiStaNodes'+str(i)+'_'+str(j)+'.Create (1);\n')
        
        		outfile.write('NetDeviceContainer staDevice'+str(i)+'_'+str(j)+' = setSta(wifiStaNodes'+str(i)+'_'+str(j)+', stack, channel, phy, wifi, mac, ssid'+str(i)+', '+str(ues[i][j][0])+', '+str(ues[i][j][1])+');\n')
        		outfile.write('clientDevices'+str(i)+'.Add(staDevice'+str(i)+'_'+str(j)+'.Get(0));\n')


		outfile.write('address.SetBase ("'+address2+'", "255.255.255.0");\n')
		outfile.write('Ipv4InterfaceContainer staInterface'+str(i)+' = address.Assign (clientDevices'+str(i)+');\n')
		outfile.write('Ipv4InterfaceContainer apInterface'+str(i)+' = address.Assign (apDevice'+str(i)+');\n')

		for j in range(len(ues[i])):
        		outfile.write('addWirelessRoutes(wifiStaNodes'+str(i)+'_'+str(j)+',apInterface'+str(i)+'.GetAddress (0));\n')
        		outfile.write('installApp(staInterface'+str(i)+'.GetAddress ('+str(j)+'), wifiStaNodes'+str(i)+'_'+str(j)+'.Get(0), wifiApNode'+str(i)+','+str(j)+','+str(i)+','+str(start)+','+str(endT)+','+str(rate)+');\n')

	end(outfile,'-wifiac')
	outfile.close()


def generate_topo_file_wifi_dynflows(outputDir,enbs,ues):
	outfile = open(outputDir + 'WiFitopo.cc','w')

	rate = 70000000
	start = 0.5
	endT = 3.5
	init(outfile,'-wifi-dyn')	# static traffic generation
	
	for i in range(len(enbs)):
		address1 = '10.1.'+str(2*i+1)+'.0'
		address2 = '10.1.'+str(2*i+2)+'.0'
		outfile.write('NodeContainer remoteHost'+str(i)+';\n')
        	outfile.write('//remoteHost'+str(i)+'.Create (1);\n')
        	outfile.write('//stack.Install (remoteHost'+str(i)+');\n')

        	outfile.write('Ssid ssid'+str(i)+' = Ssid ("ns-3-ssid'+str(i)+'");\n')

        	outfile.write('NodeContainer wifiApNode'+str(i)+';\n')
        	outfile.write('wifiApNode'+str(i)+'.Create (1);\n')
        	outfile.write('NetDeviceContainer apDevice'+str(i)+';\n')
        	outfile.write('apDevice'+str(i)+' = setAP (wifiApNode'+str(i)+', "'+address1+'", "'+address2+'", remoteHost'+str(i)+', pointToPoint, stack, address, channel, phy, wifi, mac, ssid'+str(i)+', ipv4RoutingHelper, '+str(enbs[i][0])+', '+str(enbs[i][1])+', 1);\n')

		outfile.write('NetDeviceContainer clientDevices'+str(i)+';\n')

		for j in range(len(ues[i])):
			outfile.write('NodeContainer wifiStaNodes'+str(i)+'_'+str(j)+';\n')
        		outfile.write('wifiStaNodes'+str(i)+'_'+str(j)+'.Create (1);\n')
        
        		outfile.write('NetDeviceContainer staDevice'+str(i)+'_'+str(j)+' = setSta(wifiStaNodes'+str(i)+'_'+str(j)+', stack, channel, phy, wifi, mac, ssid'+str(i)+', '+str(ues[i][j][0])+', '+str(ues[i][j][1])+');\n')
        		outfile.write('clientDevices'+str(i)+'.Add(staDevice'+str(i)+'_'+str(j)+'.Get(0));\n')


		outfile.write('address.SetBase ("'+address2+'", "255.255.255.0");\n')
		outfile.write('Ipv4InterfaceContainer staInterface'+str(i)+' = address.Assign (clientDevices'+str(i)+');\n')
		outfile.write('Ipv4InterfaceContainer apInterface'+str(i)+' = address.Assign (apDevice'+str(i)+');\n')

		for j in range(len(ues[i])):
        		outfile.write('addWirelessRoutes(wifiStaNodes'+str(i)+'_'+str(j)+',apInterface'+str(i)+'.GetAddress (0));\n')

			outfile.write('curr_flow_time = flowStartTime;// flowstartTime->GetValue ();;\n')
			outfile.write('dlport = 1;\n')
			outfile.write('while (curr_flow_time < EndFlowGeneration){\n')
			outfile.write('Webrequest (staInterface'+str(i)+'.GetAddress ('+str(j)+'), wifiStaNodes'+str(i)+'_'+str(j)+'.Get(0), wifiApNode'+str(i)+', dlport,curr_flow_time, simTime-0.1, flowId, '+str(j)+','+str(i)+', htmlObjectSize, eObjectNum, embeddedObjectSize, eObjectInterval);\n')
			outfile.write('curr_flow_time += requestInterval->GetValue ();\n')
			outfile.write('}\n')


	end(outfile,'-wifi-dyn')
	outfile.close()

def generate_topo_file_udp_dynflows(outputDir, bw,rbgsize,enbs,ues,Assignment_1,Assignment_2,Assignment_3,Assignment_4,share_1,share_2,share_3,onTime,offTime,edges,ctrlIntfrnc):
	outfile = open(outputDir + 'topo.cc','w')

	enbDevs = {}
	ueDevs = {}

		

	startTime = [0.0]
	endTime = [1.1]
	onTime = [0.1]
	offTime = [0.1]

	if (len(startTime) == 1):
		startTime = makeList(startTime[0],enbs,ues)
	if (len(endTime) == 1):
		endTime = makeList(endTime[0],enbs,ues)
	if (len(onTime) == 1):
		onTime = makeList(onTime[0],enbs,ues)
	if (len(offTime) == 1):
		offTime = makeList(offTime[0],enbs,ues)
	
	init(outfile,'-udp-dyn') # static traffic generation	
	#init(outfile,'-udp-dyn') # dynamic traffic generation	
	outfile.write('uint8_t bandwidth = '+str(bw)+'; //6:15:25:50:75:100:\n')
  	outfile.write('lteHelper->SetEnbDeviceAttribute ("DlBandwidth", UintegerValue (bandwidth));\n')
  	outfile.write('lteHelper->SetEnbDeviceAttribute ("UlBandwidth", UintegerValue (bandwidth));\n')

  	outfile.write('lteHelper->SetSchedulerType ("ns3::PfFfCellFiScheduler");\n')
  	outfile.write('lteHelper->SetFfrAlgorithmType ("ns3::LteFfrCentralizedAlgorithm");\n')

  	outfile.write('InternetStackHelper internet;\n')
  	outfile.write('Ipv4StaticRoutingHelper ipv4RoutingHelper;\n')
  	outfile.write('PointToPointHelper p2ph;\n')
  	outfile.write('Ptr<Node> remoteHost = createRemoteHost(epcHelper,internet,ipv4RoutingHelper,p2ph);\n')

  	outfile.write('uint8_t rbgSize = '+str(rbgsize)+';  // for bandwidths 1,2,2,3,4,4\n')
  	outfile.write('std::vector <bool> dl, ul, ctrl;\n')
  	outfile.write('ctrl.resize (bandwidth, false);  // false for control means, the RB is enabled\n')
  	outfile.write('dl.resize (bandwidth / rbgSize, false);  // false for dl means, the RB is enabled\n')
  	outfile.write('ul.resize (bandwidth , false); // false for ul means, the RB is enabled\n')
  	outfile.write('lteHelper->SetUlRbgMaps (ul);\n')
	#outfile.write('enum EpsBearer::Qci q = EpsBearer::GBR_CONV_VOICE;\n')
  	#outfile.write('EpsBearer bearer (q);\n')
	for i in range(len(enbs)):
		#startTime = 0
		#inc = 0
		assign_subchannels(outfile,i,Assignment_1[i],Assignment_2[i],Assignment_3[i],Assignment_4[i],share_1[i],share_2[i],share_3[i])
		gen_enb(outfile,i,enbs[i])

		enbDevs[enbs[i]] = 'enbDevs'+str(i)

		for j in range(len(ues[i])):

			ueDevs[ues[i][j]] = 'ueDevs'+str(str(i)+'_'+str(j))

			gen_ue(outfile,str(i)+'_'+str(j),ues[i][j],True)
			outfile.write('lteHelper->Attach (ueDevs'+str(str(i)+'_'+str(j))+', enbDevs'+str(i)+'.Get (0));\n')
			#outfile.write('lteHelper->ActivateDataRadioBearer (ueDevs'+str(str(i)+'_'+str(j))+', bearer);\n')

			outfile.write('curr_flow_time = flowStartTime;//flowstartTime->GetValue ();;\n')
			outfile.write('dlport = 1;\n')
			outfile.write('while (curr_flow_time < EndFlowGeneration){\n')
			outfile.write('Webrequest (ue'+str(str(i)+'_'+str(j))+', epcHelper, ipv4RoutingHelper, remoteHost, dlport,curr_flow_time, simTime-0.1, flowId, 0, 0, htmlObjectSize, eObjectNum, embeddedObjectSize, eObjectInterval,'+str(ctrlIntfrnc[ues[i][j]])+');\n')
			outfile.write('curr_flow_time += requestInterval->GetValue ();\n')
			outfile.write('}\n')



			#outfile.write('AppInstaller(lteHelper,ueDevs'+str(str(i)+'_'+str(j))+',enbDevs'+str(i)+',ue'+str(i)+'_'+str(j)+',epcHelper,ipv4RoutingHelper,remoteHost, 9,'+str(startTime[i][j])+','+str(endTime[i][j])+','+str(onTime[i][j])+','+str(offTime[i][j])+','+str(i)+','+str(j)+');\n')
			#inc += 0.1
	#print edges
	for i in range (len(enbs)):
		node = enbs[i]
		for user in ues[i]:
			outfile.write('addMyClient('+ueDevs[user]+','+ enbDevs[node]+');\n')
		for user in edges[i]:
			outfile.write('addInterferingClient('+ueDevs[user]+','+ enbDevs[node]+');\n')

	#outfile.write('fout.close()\n;')
	end(outfile)
	outfile.close()

