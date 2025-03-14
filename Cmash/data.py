from .utils import *

def retrieve_epdata(ra=None, dec=None, radius=None, start=None, end=None, url="https://ep.bao.ac.cn/ep/data_center/api/wxt_observation_data/"):
    """
    Args:
        ra (float): RA in degree. Defaults to None.
        dec (float): DEC in degree. Defaults to None.
        radius (float): radius in degree. Defaults to None.
        start (str): start time in UTC. Defaults to None.
        end (str): end time in UTC. Defaults to None.
        url (str): api url. Defaults to "https://ep.bao.ac.cn/ep/data_center/api/wxt_observation_data/".

    Raises:
        ValueError: response error

    Returns:
        _type_: response json
    """

    data = {k: v for k, v in locals().items() if v is not None}
    del data['url']

    response = requests.post(url, data=data)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError('Status code%s'%response.status_code)
    



def retrieve_gwdata(far=1/(365.25*24*3600),classification='BBH',url="https://gracedb.ligo.org/api/superevents/"):
    """
    Args:
        far (float, optional): False Alarm Rate. Defaults to 1/(365.25*24*3600).
        classification (str, optional): BBH/BNS/NSBH. Defaults to 'BBH'.
        url (str, optional): API link. Defaults to "https://gracedb.ligo.org/api/superevents/".

    Returns:
        dict: contains all event info
    """
    GRACEDB_URL = url
    # Define search parameters
    params = {
        "far": far,  # Convert 1 per year to Hz
        "classification": classification,  # Binary Black Hole classification
        "open_alert": "true",  # Only public alerts
        "format": "json",  # Return data in JSON format
        "pipeline": "O4",  # Ensure events are from the O4 observing run
    }

    event_list = dict()

    STARTS = ['?start=%s'%i for i in range(0,4400,10)]
    # Send the request to GraceDB

    for start in STARTS:
        response = requests.get(GRACEDB_URL+start, params=params)
        data = response.json()

        # Extract event information
        if "superevents" in data:
            for event in data["superevents"]:
                event_name = event["superevent_id"]
                if event_name[0] == 'M':
                    continue
                if event['far'] < far:
                    skip = 0
                    # Download the primary skymap
                    file_comb = event['links']['files']
                    link_json = requests.get(file_comb).json()
                    #Check Retraction andp_astro:
                    for linkkey,link in reversed(link_json.items()):
                        if 'Retraction' in linkkey:
                            print('%s retracted'%event_name)
                            skip = 1
                            continue
                        
                        if 'p_astro.json' in linkkey and ',' not in linkkey:
                            p_astro = requests.get(link).json()
                            if p_astro[classification] < 0.5:
                                skip = 1
                                continue
                    if skip == 0:
                        event_list[event_name] = event
                        print(f"Superevent: {event_name}, FAR: {event['far']*(365.25 * 24 * 3600)} per yr")
    return event_list