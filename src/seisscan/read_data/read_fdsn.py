from obspy import UTCDateTime, Stream, Inventory
from obspy.core import AttribDict
from obspy.clients.fdsn import Client


def read_fdsn(starttime, endtime, network, station, location, channel, provider="IRIS", attach_coordinates=True, attach_response=True):
    """
    Connects to FDSN web service of IRIS to retrive ObsPy Stream with station metadata added.
    Each Trace of the Stream object contains station coordinates and response information.
    It utilizes obspy.clients.fdsn.Client service to download waveform and inventory. For more information on the service, please visit
    "https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.html#obspy.clients.fdsn.client.Client".

    Parameters
    ----------
    starttime: obspy.UTCDateTime
        Start time for waveform and metadata.
    endtime: obspy.UTCDateTime
        End time for waveform and metadata.
    network: str
        One or more network codes seperated by comma. It supports wildcards.
    station: str
        One or more station codes seperated by comma. It supports wildcards.
    location: str
        One or more location codes seperated by comma. It supports wildcards.
    channel: str
        One or more channel codes seperated by comma. It supports wildcards.
    provider: str
        A key string for recognized FDSN server. It is one of 'IRIS', 'IRISPH5', 'GEOFON' etc.
        Please see the above link for all the providers. Default is 'IRIS'.
    attach_coordinates: bool
        If True, station coordinates are attached in each trace stats. Default is True.
    attach_response: bool)
        If True, station response information is attached to each trace stats. Default is True.
        
    Returns
    -------
    st: Obspy.Stream
        Waveform stream.
    """
    
    #--- client for FDSN web server
    client = Client(base_url=provider)
    
    #--- download obspy stream
    st = Stream()
    
    try:
        st += client.get_waveforms(network, station, location, channel, starttime, endtime)
    except:
        pass
        
    #--- download inventory
    if attach_coordinates or attach_response:
        inventory = read_fdsn_inventory(starttime, endtime, network, station, location, channel, provider=provider)
        
    #--- attach coordinates
    if attach_coordinates:
        for tr in st:
            try:
                coordinates = inventory.get_coordinates(tr.id, datetime=tr.stats.starttime)

                tr.stats.sac = AttribDict()
                tr.stats.sac.stlo = coordinates['longitude']
                tr.stats.sac.stla = coordinates['latitude']
                tr.stats.sac.stel = coordinates['elevation']
            except:
                pass
            
    #--- attach response
    if attach_response:
        for tr in st:
            try:
                tr.attach_response(inventory)
            except:
                pass
            
    return st


def read_fdsn_inventory(starttime, endtime, network, station, location, channel, provider="IRIS"):
    """
    Connects to FDSN web service of IRIS to retrive station Inventory.
    It utilizes obspy.clients.fdsn.Client service to download inventory.
    For more information on the service, please visit
    "https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.html#obspy.clients.fdsn.client.Client".
    
    Parameters
    ----------
    starttime: obspy.UTCDateTime
        Start time for waveform and metadata.
    endtime: obspy.UTCDateTime
        End time for waveform and metadata.
    network: str
        One or more network codes seperated by comma. It supports wildcards.
    station: str
        One or more station codes seperated by comma. It supports wildcards.
    location: str
        One or more location codes seperated by comma. It supports wildcards.
    channel: str
        One or more channel codes seperated by comma. It supports wildcards.
    provider: str
        A key string for recognized FDSN server. It is one of 'IRIS', 'IRISPH5', 'GEOFON' etc.
        Please see the above link for all the providers. Default is 'IRIS'.
        
    Returns
    -------
    inventory: Obspy.Inventory
    """
    
    #--- client for FDSN web server
    client = Client(base_url=provider)
    
    #--- download inventory
    inventory = Inventory()
    
    try:
        inventory += client.get_stations(starttime=starttime, endtime=endtime,
                                        network=network, station=station, location=location, channel=channel,
                                        level='response')
    except:
        pass
        
    return inventory