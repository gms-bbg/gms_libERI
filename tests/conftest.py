def pytest_addoption(parser):
    parser.addoption("--plugindir", action="store", default="../build")

def pytest_generate_tests(metafunc):
    option_value = metafunc.config.option.plugindir
    if 'plugindir' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("plugindir", [option_value])
