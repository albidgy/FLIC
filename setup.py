from setuptools import setup


setup(
    name='FLIC',
    version='1.0',
    author='Alexandra Kasianova',
    author_email='Alexandra.Kasianova@skoltech.ru',
    url='https://github.com/albidgy/FLIC',
    description='A tool for isoform reconstruction based on long reads',
    license='MIT',
    install_requires=['joblib', 'networkx', 'numpy'],
    packages=['flic_src', 'flic_src/modules', 'flic_src/scripts'],
    scripts=['flic_src/modules/run_cagefightr.R'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'flic=flic_src.run_flic:run_tool',
        ]
    },
)
