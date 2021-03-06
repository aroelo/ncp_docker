from flask.cli import FlaskGroup
from project import app, db
from project.download_pavian_data.models import User
from project.models import PavianInput
import os
from datetime import datetime

cli = FlaskGroup(app)


@cli.command("drop_db")
def drop_db():
    db.drop_all()
    db.session.commit()


@cli.command("create_db")
def create_db():
    db.create_all()
    db.session.commit()


@cli.command("seed_pavian_db")
def seed_pavian_db():
    input_dir = app.config['PAVIAN_IN']
    # pavian files in sql
    pavian_sql = pavian_files_sql()
    # pavian files in dir
    pavian_dir, pavian_symlink_files = pavian_files_dir(input_dir)

    # pavian files that are not yet in sql, but in directory (should be added to sql)
    not_in_sql = pavian_dir - pavian_sql
    # pavian files that are in sql, but not in directory (should be deleted from sql)
    not_in_dir = pavian_sql - pavian_dir

    print("Added new runs:")
    for pavian_basename in not_in_sql:
        print(pavian_basename.strip('.pavian'))
        insert_into_sql(input_dir, pavian_basename)

    print("\nUpdated 'support' column for runs. No support files for:")
    for pavian_basename in pavian_dir:
        project_name = pavian_basename.strip('.pavian')
        if not os.path.exists(os.path.join(input_dir, 'support_files', project_name)):
            print(project_name)
            update_support_column(pavian_basename, False)
        else:
            update_support_column(pavian_basename, True)

    print("\nIgnored pavian symlink files:")
    for pavian_basename in pavian_symlink_files:
        print(pavian_basename.strip('.pavian'))

    print("\nDeleted runs:")
    for pavian_basename in not_in_dir:
        print(pavian_basename.strip('.pavian'))
        delete_entry(pavian_basename)


def pavian_files_sql():
    """
    List pavian files in sql database
    :return: set of list of pavian files
    """
    pavian_files = [record.file for record in PavianInput.query.all()]
    return set(pavian_files)


def pavian_files_dir(input_dir):
    """
    List pavian files in input directory
    :param input_dir: directory containing pavian files
    :return: set of list of pavian files
    """
    (root, dirs, files) = next(os.walk(input_dir))
    pavian_files = []
    pavian_symlink_files = []
    for file in files:
        if file.endswith('.pavian'):
            if os.path.islink(os.path.join(root, file)):
                pavian_symlink_files.append(file)
            else:
                pavian_files.append(file)
    return set(pavian_files), set(pavian_symlink_files)


def insert_into_sql(input_dir, pavian_basename):
    """
    Insert data into mysql database
    :param input_dir: directory containing pavian files
    :param pavian_basename: relative pavian file name
    :return: None
    """
    if len(pavian_basename.split('.')) < 4:
        run = pavian_basename.split('.')[0]
        sample = ''
        nt = ''
    else:
        run, sample, nt, _ = pavian_basename.split('.')
    date = pavian_basename.split('_')[0]
    try:
        datetime.strptime(str(date), '%Y%m%d')
    except ValueError:
        # date = ''
        date = '20010101'
    support = False
    pavian_file = open(os.path.join(input_dir, pavian_basename), "r")
    version = None
    for line in pavian_file:
        line = line.strip('\n')
        if line == 'Naktuinbouw':
            version = 'Naktuinbouw'
            continue
        if version == 'Naktuinbouw':
            organism_taxid = line.split('\t')[6]
            organism_name = line.split('\t')[7].strip(' ')
        else:
            organism_taxid = line.split('\t')[4]
            organism_name = line.split('\t')[5].strip(' ')
        pavian_row = PavianInput(pavian_basename, run, sample, nt, date, organism_taxid, organism_name, support)
        db.session.add(pavian_row)
    db.session.commit()


def update_support_column(pavian_basename, column_value):
    """
    Update column that describes whether 'support files' are present in the 'support_files' directory.
    :param pavian_basename: relative pavian file name
    :param column_value:  column value, 0 if false/no support files, 1 if true
    :return: None
    """
    pavian_row = PavianInput.query.filter_by(file=pavian_basename).update(dict(support=column_value))
    # pavian_row.support = column_value
    db.session.commit()


def delete_entry(pavian_basename):
    """
    Delete an entry from the mysql database if the pavian file was once available, but is now not in the directory anymore.
    :param pavian_basename: relative pavian file name
    :return: None
    """
    PavianInput.query.filter_by(file=pavian_basename).delete()
    db.session.commit()


if __name__ == "__main__":
    cli()
