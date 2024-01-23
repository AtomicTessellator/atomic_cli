import uuid


def is_valid_uuid(uuid_string):
    try:
        uuid_obj = uuid.UUID(uuid_string, version=4)
        return str(uuid_obj) == uuid_string
    except ValueError:
        return False
