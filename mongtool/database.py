import pymongo

class Database(object):
    uri = "mongodb://localhost:27017/"
    db = None

    @staticmethod
    def initialize(db_name):
        client = pymongo.MongoClient(Database.uri)
        Database.db = client[db_name] # Database Name
        Database.db_name = db_name # Database Name
        #Database.collection = Database.db["sample"] # Collection Name

    @staticmethod
    def insert(collection, data):
        Database.db[collection].insert(data)

    @staticmethod
    def find(collection, query):
        return Database.db[collection].find(query)

    @staticmethod
    def find_one(collection, query):
        return Database.db[collection].find_one(query)
    
    @staticmethod
    def get_pvl(collection, query):
        return Database.db[collection].find(query, {"_id": 0, "aribavir.lukS_PV.present": 1})
    
    @staticmethod
    def get_mlst(collection, query):
        return Database.db[collection].find(query, {"_id": 0, "mlst": 1})
    
    @staticmethod
    def get_cgmlst(collection, query):
        return Database.db[collection].find(query, {"_id": 0, "alleles": 1})
