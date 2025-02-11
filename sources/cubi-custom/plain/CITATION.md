# Citing this repository

If you are using the content of this repository in whole or in part for your own work,
please credit the Core Unit Bioinformatics in an appropriate form.

In general, please add this statement to the acknowledgments of your publication:

    This work was supported by the Core Unit Bioinformatics
    of the Medical Faculty of the Heinrich Heine University Düsseldorf.

Additionally, please follow the below instructions to obtain a citable reference
for your publication.

## Identifying the right source link

Each repository of the Core Unit Bioinformatics is assigned a persistent
identifier (PID) at some point (usually after the prototype stage). Please use
this PID to link to the repository. You always find the PID in the top-level
`pyproject.toml`. Depending on the type of repository (project, workflow, or
workflow template), the relevant PID is listed in the corresponding metadata
section:

```toml
# workflow repository
[cubi.workflow]
pid = "THE-PID"
```

```toml
# workflow template repository
[cubi.workflow.template]
pid = "THE-PID"
```

```toml
# project repository
[cubi.project]
pid = "THE-PID"
```

If a PID has not yet been assigned to the repository, please use the repository URL,
and, if time permits, contact the repository maintainer regarding assigning a PID
in the near future.

### 1. Release version

For release versions, please use the respective version string in addition to the source link,
i.e. ideally in combination with the PID, and integrate that information into your list
of references as appropriate.

Note that repositories of the type "project" may not contain a lot of code, and are thus
often not amenable to the usual "release cycle" following bug fixes, feature integrations and so on.
Hence, the "project" metadata do not contain a "version" key (as opposed to workflow and workflow
template repositories). See the next point if you encounter that situation.

### 2. Development (non-release) version

For development versions, please use the respective git commit hash in addtion to the source link,
i.e. ideally in combination with the PID, and integrate that information into your list
of references as appropriate. It is strongly recommended to only use git commits from the two
central branches `main` and `dev`.

If a "project" repository is lacking an explicit release version, please use the same strategy
to obtain a citable reference of the repository.


### 3. None of the above

Please get in touch and we'll find a solution for your case :-)
