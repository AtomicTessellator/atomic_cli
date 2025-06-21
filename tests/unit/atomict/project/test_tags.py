from unittest.mock import MagicMock, patch

import pytest

from atomict.project.tags import (
    VALID_TAG_COLOURS,
    create_project_tag,
    create_tag,
    delete_project_tag,
    delete_project_tag_association,
    get_project_tag,
    get_tag_by_name,
    list_project_tag_associations,
    list_project_tags,
    project_tag_exists,
    tag_exists,
    update_project_tag,
)


class TestCreateTag:
    @patch("atomict.project.tags.post")
    def test_create_tag_valid_color(self, mock_post):
        mock_post.return_value = {"id": 1}

        result = create_tag("test-tag", "bg-success")

        mock_post.assert_called_once_with(
            "api/project-tag/",
            {"tag": "test-tag", "color": "bg-success"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": 1}

    def test_create_tag_invalid_color(self):
        with pytest.raises(ValueError) as exc_info:
            create_tag("test-tag", "invalid-color")

        assert "Invalid tag color: invalid-color" in str(exc_info.value)
        assert "choose from" in str(exc_info.value)

    @patch("atomict.project.tags.post")
    def test_create_tag_all_valid_colors(self, mock_post):
        mock_post.return_value = {"id": 1}

        for color in VALID_TAG_COLOURS:
            create_tag("test-tag", color)
            mock_post.assert_called_with(
                "api/project-tag/",
                {"tag": "test-tag", "color": color},
                extra_headers={"Content-Type": "application/json"},
            )


class TestGetTagByName:
    @patch("atomict.project.tags.get")
    def test_get_tag_by_name_success(self, mock_get):
        mock_get.return_value = {
            "results": [{"id": 1, "tag": "test-tag", "color": "bg-success"}]
        }

        result = get_tag_by_name("test-tag")

        mock_get.assert_called_once_with("api/project-tag/?tag=test-tag")
        assert result == {"id": 1, "tag": "test-tag", "color": "bg-success"}


class TestTagExists:
    @patch("atomict.project.tags.get")
    def test_tag_exists_true(self, mock_get):
        mock_get.return_value = {"count": 1}

        result = tag_exists("existing-tag")

        mock_get.assert_called_once_with("api/project-tag/?tag=existing-tag")
        assert result is True

    @patch("atomict.project.tags.get")
    def test_tag_exists_false(self, mock_get):
        mock_get.return_value = {"count": 0}

        result = tag_exists("non-existing-tag")

        mock_get.assert_called_once_with("api/project-tag/?tag=non-existing-tag")
        assert result is False


class TestCreateProjectTag:
    @patch("atomict.project.tags.post")
    def test_create_project_tag_success(self, mock_post):
        mock_post.return_value = {"id": 1}

        result = create_project_tag("project-123", "tag-456")

        mock_post.assert_called_once_with(
            "api/project-tag-project/",
            {"project": "project-123", "project_tag": "tag-456"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": 1}


class TestProjectTagExists:
    @patch("atomict.project.tags.get")
    def test_project_tag_exists_true(self, mock_get):
        mock_get.return_value = {"count": 1}

        result = project_tag_exists("project-123", "tag-456")

        mock_get.assert_called_once_with(
            "api/project-tag-project/?project=project-123&project_tag=tag-456"
        )
        assert result is True

    @patch("atomict.project.tags.get")
    def test_project_tag_exists_false(self, mock_get):
        mock_get.return_value = {"count": 0}

        result = project_tag_exists("project-123", "tag-456")

        mock_get.assert_called_once_with(
            "api/project-tag-project/?project=project-123&project_tag=tag-456"
        )
        assert result is False


class TestDeleteProjectTag:
    @patch("atomict.project.tags.delete")
    def test_delete_project_tag_success(self, mock_delete):
        mock_delete.return_value = {}

        result = delete_project_tag("tag-123")

        mock_delete.assert_called_once_with("api/project-tag/tag-123/")
        assert result == {}


class TestListProjectTags:
    @patch("atomict.project.tags.get")
    def test_list_project_tags_success(self, mock_get):
        mock_get.return_value = {
            "results": [
                {"id": 1, "tag": "tag1", "color": "bg-success"},
                {"id": 2, "tag": "tag2", "color": "bg-danger"},
            ]
        }

        result = list_project_tags()

        mock_get.assert_called_once_with("api/project-tag/")
        assert len(result["results"]) == 2


class TestGetProjectTag:
    @patch("atomict.project.tags.get")
    def test_get_project_tag_success(self, mock_get):
        mock_get.return_value = {"id": 1, "tag": "test-tag", "color": "bg-success"}

        result = get_project_tag("tag-123")

        mock_get.assert_called_once_with("api/project-tag/tag-123/")
        assert result == {"id": 1, "tag": "test-tag", "color": "bg-success"}


class TestUpdateProjectTag:
    @patch("atomict.project.tags.patch")
    def test_update_project_tag_name_only(self, mock_patch):
        mock_patch.return_value = {"id": 1, "tag": "new-name", "color": "bg-success"}

        result = update_project_tag("tag-123", tag="new-name")

        mock_patch.assert_called_once_with(
            "api/project-tag/tag-123/",
            {"tag": "new-name"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": 1, "tag": "new-name", "color": "bg-success"}

    @patch("atomict.project.tags.patch")
    def test_update_project_tag_color_only(self, mock_patch):
        mock_patch.return_value = {"id": 1, "tag": "test-tag", "color": "bg-danger"}

        result = update_project_tag("tag-123", color="bg-danger")

        mock_patch.assert_called_once_with(
            "api/project-tag/tag-123/",
            {"color": "bg-danger"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": 1, "tag": "test-tag", "color": "bg-danger"}

    @patch("atomict.project.tags.patch")
    def test_update_project_tag_both_parameters(self, mock_patch):
        mock_patch.return_value = {"id": 1, "tag": "new-name", "color": "bg-warning"}

        result = update_project_tag("tag-123", tag="new-name", color="bg-warning")

        mock_patch.assert_called_once_with(
            "api/project-tag/tag-123/",
            {"tag": "new-name", "color": "bg-warning"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": 1, "tag": "new-name", "color": "bg-warning"}

    def test_update_project_tag_invalid_color(self):
        with pytest.raises(ValueError) as exc_info:
            update_project_tag("tag-123", color="invalid-color")

        assert "Invalid tag color: invalid-color" in str(exc_info.value)
        assert "choose from" in str(exc_info.value)

    def test_update_project_tag_no_parameters(self):
        with pytest.raises(ValueError) as exc_info:
            update_project_tag("tag-123")

        assert "At least one parameter (tag or color) must be provided" in str(
            exc_info.value
        )

    @patch("atomict.project.tags.patch")
    def test_update_project_tag_all_valid_colors(self, mock_patch):
        mock_patch.return_value = {"id": 1, "tag": "test-tag", "color": "bg-success"}

        for color in VALID_TAG_COLOURS:
            update_project_tag("tag-123", color=color)
            mock_patch.assert_called_with(
                "api/project-tag/tag-123/",
                {"color": color},
                extra_headers={"Content-Type": "application/json"},
            )


class TestDeleteProjectTagAssociation:
    @patch("atomict.project.tags.delete")
    def test_delete_project_tag_association_success(self, mock_delete):
        mock_delete.return_value = {}

        result = delete_project_tag_association("association-123")

        mock_delete.assert_called_once_with("api/project-tag-project/association-123/")
        assert result == {}


class TestListProjectTagAssociations:
    @patch("atomict.project.tags.get")
    def test_list_project_tag_associations_no_filters(self, mock_get):
        mock_get.return_value = {"results": []}

        result = list_project_tag_associations()

        mock_get.assert_called_once_with("api/project-tag-project/")
        assert result == {"results": []}

    @patch("atomict.project.tags.get")
    def test_list_project_tag_associations_project_filter(self, mock_get):
        mock_get.return_value = {"results": [{"id": 1}]}

        result = list_project_tag_associations(project_id="project-123")

        mock_get.assert_called_once_with("api/project-tag-project/?project=project-123")
        assert result == {"results": [{"id": 1}]}

    @patch("atomict.project.tags.get")
    def test_list_project_tag_associations_tag_filter(self, mock_get):
        mock_get.return_value = {"results": [{"id": 2}]}

        result = list_project_tag_associations(tag_id="tag-456")

        mock_get.assert_called_once_with("api/project-tag-project/?project_tag=tag-456")
        assert result == {"results": [{"id": 2}]}

    @patch("atomict.project.tags.get")
    def test_list_project_tag_associations_both_filters(self, mock_get):
        mock_get.return_value = {"results": [{"id": 3}]}

        result = list_project_tag_associations(
            project_id="project-123", tag_id="tag-456"
        )

        mock_get.assert_called_once_with(
            "api/project-tag-project/?project=project-123&project_tag=tag-456"
        )
        assert result == {"results": [{"id": 3}]}


class TestValidTagColours:
    def test_valid_tag_colours_constants(self):
        """Test that VALID_TAG_COLOURS contains expected Bootstrap classes"""
        expected_colors = [
            "bg-success",
            "bg-primary",
            "bg-secondary",
            "bg-danger",
            "bg-warning",
            "bg-info",
            "bg-light",
            "bg-dark",
        ]

        assert VALID_TAG_COLOURS == expected_colors
        assert len(VALID_TAG_COLOURS) == 8
