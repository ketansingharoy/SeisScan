import os
import pandas as pd
from obspy import UTCDateTime, Stream, read_inventory, Inventory
from obspy.core import AttribDict
from obspy.clients.fdsn.mass_downloader import RectangularDomain, Restrictions, MassDownloader
from obspy.clients.filesystem.sds import Client


def get_station_coordinates_from_inventory_(inv, datetime):
    
    ################### Prepare station file ##############
    #--- loop over seed_id to extract information
    network_list = []
    station_list = []
    longitude_list = []
    latitude_list = []
    elevation_list = []
    
    _, _, seed_id_list = inv.get_contents().values()

    for seed_id in seed_id_list:

        seed_id = seed_id.strip()
        network, station, _, _ = seed_id.split('.')

        latitude, longitude, elevation, _ = inv.get_coordinates(seed_id, datetime=datetime).values()

        network_list.append(network)
        station_list.append(station)
        longitude_list.append(longitude)
        latitude_list.append(latitude)
        elevation_list.append(elevation)


    #--- Prepare dictionary containing station information
    station_dict = {'network': network_list, 'station': station_list, 'longitude':longitude_list, 'latitude':latitude_list, 'elevation':elevation_list}

    #--- Prepare station dataframe and write csv
    df = pd.DataFrame(station_dict)
    df = df.drop_duplicates('station').reset_index(drop=True)
    
    return df

def make_sds_from_fdsn(provider, starttime, endtime, network, station, location, channel,
                minlongitude, maxlongitude, minlatitude, maxlatitude, db_dir,
                chunklength_in_sec=86400.0, reject_channels_with_gaps=False, minimum_length=0.0,
                minimum_interstation_distance_in_m=0.0):
    
    '''
    Download data from FDSN peovider.
    Prepares waveform archive.
    Prepares stationxml archive.
    '''
    
    ##########################
    def get_mseed_storage(network, station, location, channel, starttime, endtime):
        '''
        Function to check if a mseed file exists.
        Required by prepare_sds function.
        '''
        #---- The type of data is waveform
        data_type = 'D'
        
        #--- get year and julian day from starttime
        year_str = starttime.strftime('%Y')
        jday_str = starttime.strftime('%j')

        #--- prepare path of mseed file
        wfile_ = f'{network}.{station}.{location}.{channel}.{data_type}.{year_str}.{jday_str}'
        wfile = os.path.join(waveform_dir, year_str, network, station, f'{channel}.{data_type}', wfile_)

        #--- Check if the mseed file exists
        if os.path.exists(wfile):
            return True
        else:
            return wfile
    ##########################
    
    #--- check db_dir
    if not os.path.exists(db_dir):
        os.mkdir(db_dir)
        
    #---- create directory to store mseed files
    waveform_dir = os.path.join(db_dir, 'waveform')
    
    if not os.path.exists(waveform_dir):
        os.mkdir(waveform_dir)
        
    #---- create directory to store stationxml files
    metadata_dir = os.path.join(db_dir, 'metadata')
    
    if not os.path.exists(metadata_dir):
        os.mkdir(metadata_dir)
        
    #--- form dimain of rectangular search
    domain = RectangularDomain(minlatitude=minlatitude, maxlatitude=maxlatitude, minlongitude=minlongitude, maxlongitude=maxlongitude)

    #--- form restriction
    restrictions = Restrictions(
        starttime, endtime,
        network=network,
        station=station,
        location=location,
        channel=channel,
        chunklength_in_sec=chunklength_in_sec,
        reject_channels_with_gaps=reject_channels_with_gaps,
        minimum_length=minimum_length,
        minimum_interstation_distance_in_m=minimum_interstation_distance_in_m
    )

    #--- prepare mass downloader and start downloading
    mdl = MassDownloader(providers=[provider])
    mdl.download(domain, restrictions, mseed_storage=get_mseed_storage, stationxml_storage=metadata_dir)
    
    #--- prepare station file
    mfile = os.path.join(metadata_dir, '*.xml')
    inv = read_inventory(mfile)
    
    station_file = os.path.join(db_dir, 'station.txt')
    df_station = get_station_coordinates_from_inventory_(inv, starttime)
    df_station.to_csv(station_file, index=None)
    
def get_inv_from_db(db_dir, network, stationString, location, channel, starttime, endtime):
    
    #--- stationxml directories
    metadata_dir = os.path.join(db_dir, 'metadata')
    
    #--- get stations
    stations = stationString.split(',')
    
    #--- loop over station
    inv = Inventory()
    
    for station in stations:
        mfile = os.path.join(metadata_dir, f'{network}.{station}.xml')
        inventory = read_inventory(mfile)
        inv += inventory.select(network=network, station=station, location=location, channel=channel, starttime=starttime, endtime=endtime)
        
    return inv
    
def read_sds(db_dir, network, stationString, location, channel, starttime, endtime, merge=-1):

    #--- waveform and stationxml directories
    waveform_dir = os.path.join(db_dir, 'waveform')
    metadata_dir = os.path.join(db_dir, 'metadata')
    
    #--- get stations
    stations = stationString.split(',')

    #--- sds client
    client = Client(waveform_dir)

    #--- read stream: loop over station
    st = Stream()
    
    for station in stations:
        st += client.get_waveforms(network, station, location, channel, starttime, endtime, merge=merge)

    #--- get inventory
    inv = get_inv_from_db(db_dir, network, stationString, location, channel, starttime, endtime)
    
    #--- loop over trace to input coordinate and response information
    for tr in st:

        #--- Get coordinates and put them in the trace
        stla, stlo, stel, _ = inv.get_coordinates(tr.id, datetime=tr.stats.starttime).values()

        tr.stats.sac = AttribDict()
        tr.stats.sac.stlo = stlo
        tr.stats.sac.stla = stla
        tr.stats.sac.stel = stel
        
        tr.attach_response(inv)
        
    return st
