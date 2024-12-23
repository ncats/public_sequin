pipeline {
    options {
        timestamps()
        skipDefaultCheckout()
    }
    agent {
        node { label 'sctl-dev-ec2-01'}
    }
    parameters {
        string(name: 'BUILD_VERSION', defaultValue: '', description: 'The build version to deploy (optional)')
    }
    triggers {
        pollSCM('H/5 * * * *')
    }
    environment {
        PROJECT_NAME = "sctl-rshiny-complex"
    }
    stages {
        stage('Build Version') {
            when {
                expression {
                    return !params.BUILD_VERSION
                }
            }
            steps{
                script {
                    BUILD_VERSION_GENERATED = VersionNumber(
                        versionNumberString: 'v${BUILD_YEAR, XX}.${BUILD_MONTH, XX}${BUILD_DAY, XX}.${BUILDS_TODAY}',
                        projectStartDate:    '1970-01-01',
                        skipFailedBuilds:    true)
                    currentBuild.displayName = BUILD_VERSION_GENERATED
                    env.BUILD_VERSION = BUILD_VERSION_GENERATED
                    env.BUILD = 'true'
                }
            }
        }
        stage('BuildS') {
            parallel {
                stage('AMD') {
                    when {
                        expression {
                            return !params.BUILD_VERSION
                    }
                }
                agent {
                    node { label 'sctl-dev-ec2-01'}
                }
                steps {
                    sshagent (credentials: ['871f96b5-9d34-449d-b6c3-3a04bbd4c0e4']) {
                        withEnv([
                            "IMAGE_NAME=sctl-rshiny-complex",
                                "BUILD_VERSION=" + (params.BUILD_VERSION ?: env.BUILD_VERSION)
                        ]) {
                            checkout scm
                                script {
                                    docker.withRegistry("https://registry-1.docker.io/v2/","f16c74f9-0a60-4882-b6fd-bec3b0136b84") {
                                        def image = docker.build(
                                            "ncats/sctl-rshiny-complex:${env.BUILD_VERSION}",
                                            "--no-cache ."
                                        )
                                        image.push("${env.BUILD_VERSION}")
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        stage('Deploy Application') {
            agent {
                node { label 'sctl-dev-ec2-01'}
            }
            steps {
                cleanWs()
                configFileProvider([
                    configFile(fileId: 'sql.R', targetLocation: 'sql.R'),
                    configFile(fileId: 'sctl-docker-compose.yml', targetLocation: 'docker-compose.yml')
                ]) {
                    withEnv([
                        "DOCKER_REPO_NAME=ncats/sctl-rshiny-complex",
                        "BUILD_VERSION=" + (params.BUILD_VERSION ?: env.BUILD_VERSION)
                    ]) {
                        sh 'chmod 755 sql.R'
                        script {
                            def docker = new org.labshare.Docker()
                            docker.deployDockerUI()
                        }
                    }
                }
            }
        }
    }
}
